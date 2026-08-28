!> @file atlas_fesom_mesh.F90
!! @brief Atlas-backed mesh setup wrapper for tracer dwarf
!! @details When ENABLE_ATLAS is enabled, this module initializes Atlas
!!          and generates a mesh-driven configuration path for the tracer
!!          dwarf. The resulting FESOM mesh arrays are populated through the
!!          existing in-memory analytic mesh generator. When Atlas is disabled,
!!          it falls back to the standard mesh_setup() file-based path.

module atlas_fesom_mesh_module
  use MOD_MESH
  use MOD_PARTIT
  use o_PARAM, only: WP, MP
  use oce_mesh_module, only: mesh_setup, test_tri, find_levels_min_e2n, &
                              mesh_areas, mesh_auxiliary_arrays
  use analytic_mesh_module, only: generate_analytic_mesh
  use g_config, only: force_rotation, MeshPath
  use g_rotate_grid, only: set_mesh_transform_matrix, g2r
  use par_support_interfaces, only: init_mpi_types, init_gatherLists
  use iso_fortran_env, only: output_unit

  use mpi
#ifdef ENABLE_ATLAS
  use atlas_module
  use, intrinsic :: iso_c_binding, only: c_int
#endif

  implicit none
  private

#if defined(USE_HALF_PRECISION) || defined(USE_SINGLE_PRECISION)
  integer, parameter :: MPI_MP = MPI_REAL
#else
  integer, parameter :: MPI_MP = MPI_DOUBLE_PRECISION
#endif

  public :: mesh_setup_with_atlas, compute_tracer_stats_atlas
  public :: atlas_fesom_enabled, atlas_fesom_active
  public :: atlas_halo_exchange_nodal
#ifdef ENABLE_ATLAS
  public :: atlas_mesh_to_fesom_mesh, compute_field_stats_atlas
  ! Module-level Atlas mesh (persists for stats computation)
  type(atlas_Mesh), save :: atlas_mesh_global
  logical, save :: atlas_mesh_active = .false.
#endif

contains

  logical function atlas_fesom_enabled()
    character(len=32) :: value
    integer :: status
    logical, save :: checked = .false.
    logical, save :: enabled = .false.

#ifdef ENABLE_ATLAS
    if (.not. checked) then
      call get_environment_variable('ATLAS_FESOM', value, status=status)
      enabled = status == 0 .and. trim(value) == '1'
      checked = .true.
    end if
    atlas_fesom_enabled = enabled
#else
    atlas_fesom_enabled = .false.
#endif
  end function atlas_fesom_enabled

  logical function atlas_fesom_active()
#ifdef ENABLE_ATLAS
    atlas_fesom_active = atlas_mesh_active
#else
    atlas_fesom_active = .false.
#endif
  end function atlas_fesom_active

  subroutine atlas_halo_exchange_nodal(nodal_data)
    real(kind=WP), intent(inout) :: nodal_data(:,:)
#ifdef ENABLE_ATLAS
    type(atlas_functionspace_NodeColumns) :: fs
    type(atlas_Field) :: field
    real(WP), pointer :: field_data(:,:)

    ! todo: centralise fs
    fs = atlas_functionspace_NodeColumns(atlas_mesh_global, halo=2)

    field = atlas_Field('tmp',nodal_data)
    call fs%halo_exchange(field)

    call field%final()
    call fs%final()
#endif
  end subroutine atlas_halo_exchange_nodal

  !> @brief Mesh setup entry point with optional Atlas path
  !! @param[inout] partit Partition structure
  !! @param[inout] mesh Mesh structure
  subroutine mesh_setup_with_atlas(partit, mesh)
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(inout), target :: mesh

#ifdef ENABLE_ATLAS
    type(atlas_Grid)             :: grid_obj
    type(atlas_MeshGenerator)   :: meshgen_obj
    type(atlas_Mesh)            :: mesh_obj
    type(atlas_GridDistribution) :: dist_obj

    integer :: io_stat
    character(len=32)  :: use_fesom_dist_str
    character(len=10)  :: npes_string
    character(len=256) :: dist_mesh_dir, rpart_file, owner_file
    character(len=5) :: rank_string
    integer :: ncells, nnodes
    integer :: nl_default, ncells_shrinked
    real(kind=MP), allocatable :: zbar_default(:)
    type(atlas_mesh_Nodes) :: mesh_nodes_obj
    type(atlas_mesh_Cells) :: mesh_cells_obj
    type(atlas_functionspace_NodeColumns) :: fs_dummy
    ! Partition file variables
    integer :: funit, r, node_idx, npes_in_file, error_status, ierror
    integer :: file_rank, owned_count, ghost_count
    integer, allocatable :: part_csr(:), part_per_node(:), local_node_list(:)
    logical :: have_dist

    if (.not. atlas_fesom_enabled()) then
      if (partit%mype == 0) then
        write(output_unit, '(A)') &
          '  --> ATLAS_FESOM is not 1; using standard FESOM mesh setup'
      end if
      call mesh_setup(partit, mesh)
      return
    end if

    ! Use sensible defaults for fesom-pi grid
    nl_default = 10

    ! Read the vertical profile from aux3d.out, matching read_mesh.
    if (partit%mype == 0) then
      open(newunit=funit, file=trim(MeshPath)//'aux3d.out', status='old', &
           action='read', iostat=io_stat)
      if (io_stat == 0) then
        read(funit, *, iostat=io_stat) nl_default
        if (io_stat == 0 .and. nl_default >= 3) then
          allocate(zbar_default(nl_default))
          read(funit, *, iostat=io_stat) zbar_default
        end if
        close(funit)
      end if
      if (io_stat /= 0 .or. nl_default < 3) then
        nl_default = 10
        if (allocated(zbar_default)) deallocate(zbar_default)
        allocate(zbar_default(nl_default))
        zbar_default = [( -1000.0_MP * real(r-1, MP) / real(nl_default-1, MP), &
                          r=1,nl_default )]
      end if
    end if
    if (partit%npes > 1) then
      call MPI_BCast(nl_default, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      if (partit%mype /= 0) allocate(zbar_default(nl_default))
      call MPI_BCast(zbar_default, nl_default, MPI_MP, 0, &
                     partit%MPI_COMM_FESOM, ierror)
    end if
    if (zbar_default(2) > 0.0_MP) zbar_default = -zbar_default

    ! Initialize Atlas (loads plugins including atlas-fesom, registers grids)
    call atlas_initialize()

    ! Attempt to create fesom-pi grid; fall back to mesh_setup if unavailable
    if (partit%mype == 0) then
      write(output_unit, '(A)') '  --> Attempting to use Atlas fesom-pi grid...'
    end if

    grid_obj = atlas_Grid("fesom-pi")
    write(output_unit, *) '  --> Atlas fesom-pi grid found: ', grid_obj%name()

    ! Optionally drive Atlas mesh distribution from the standard FESOM owner lists.
    ! Set ATLAS_USE_FESOM_DIST=1 to enable; rpart.out and my_list*.out are read
    ! from the same dist_<npes>/ directory as the non-Atlas mesh_setup() path.
    have_dist = .false.
    call get_environment_variable('ATLAS_USE_FESOM_DIST', use_fesom_dist_str, &
                                  status=io_stat)
    if (io_stat == 0 .and. trim(use_fesom_dist_str) == '1') then
      write(npes_string, "(I10)") partit%npes
      dist_mesh_dir = trim(MeshPath)//'dist_'//trim(ADJUSTL(npes_string))//'/'
      rpart_file    = trim(dist_mesh_dir)//'rpart.out'
      error_status  = 0
      if (partit%mype == 0) then
        open(newunit=funit, file=trim(rpart_file), action='read', status='old', &
             iostat=io_stat)
        if (io_stat /= 0) then
          write(output_unit, '(3A)') &
            '  WARNING: ATLAS_USE_FESOM_DIST=1 but cannot open ', &
            trim(rpart_file), ' — using default Atlas distribution'
          error_status = 1
        else
          allocate(part_csr(partit%npes+1))
          read(funit, *) npes_in_file
          if (npes_in_file /= partit%npes) error_status = 1
          part_csr(1) = 1
          read(funit, *) part_csr(2:partit%npes+1)
          ! Convert per-partition counts to cumulative offsets (same as read_mesh)
          do r = 2, partit%npes+1
            part_csr(r) = part_csr(r-1) + part_csr(r)
          end do
          close(funit)
        end if
      end if
      if (partit%npes > 1) then
        call MPI_BCast(error_status, 1, MPI_INTEGER, 0, &
                       partit%MPI_COMM_FESOM, ierror)
        if (error_status == 0) then
          if (partit%mype /= 0) allocate(part_csr(partit%npes+1))
          call MPI_BCast(part_csr, partit%npes+1, MPI_INTEGER, 0, &
                         partit%MPI_COMM_FESOM, ierror)
        end if
      end if
      if (error_status == 0) then
        allocate(part_per_node(int(grid_obj%size())))
        part_per_node = -1
        if (partit%mype == 0) then
          do r = 0, partit%npes - 1
            write(rank_string, '(I5.5)') r
            owner_file = trim(dist_mesh_dir)//'my_list'//rank_string//'.out'
            open(newunit=funit, file=trim(owner_file), action='read', status='old', &
                 iostat=io_stat)
            if (io_stat /= 0) error stop 'Cannot open FESOM partition owner list'
            read(funit, *) file_rank
            read(funit, *) owned_count
            read(funit, *) ghost_count
            if (file_rank /= r) error stop 'Unexpected rank in FESOM partition owner list'
            allocate(local_node_list(owned_count + ghost_count))
            read(funit, *) local_node_list
            close(funit)
            do node_idx = 1, owned_count
              part_per_node(local_node_list(node_idx)) = r
            end do
            deallocate(local_node_list)
          end do
          if (any(part_per_node < 0)) error stop 'Incomplete FESOM node ownership map'
        end if
        call MPI_BCast(part_per_node, int(grid_obj%size()), MPI_INTEGER, 0, &
                       partit%MPI_COMM_FESOM, ierror)
        deallocate(part_csr)
        dist_obj   = atlas_GridDistribution(part_per_node)
        deallocate(part_per_node)
        have_dist  = .true.
        if (partit%mype == 0) write(output_unit, '(3A,I0,A)') &
          '  --> ATLAS_USE_FESOM_DIST: using partition from ', &
          trim(rpart_file), ' (', partit%npes, ' partitions)'
      else
        if (allocated(part_csr)) deallocate(part_csr)
      end if
    end if

    ! Create mesh generator and generate mesh
    meshgen_obj = atlas_MeshGenerator("fesom")
    if (have_dist) then
      mesh_obj = meshgen_obj%GENERATE(grid_obj, dist_obj)
      call dist_obj%final()
    else
      mesh_obj = meshgen_obj%GENERATE(grid_obj)
    end if

    ! This adapts the mesh with one halo
    fs_dummy = atlas_functionspace_NodeColumns(mesh_obj,halo=2)

    ! Get mesh dimensions for logging
    mesh_nodes_obj = mesh_obj%nodes()
    mesh_cells_obj = mesh_obj%cells()
    nnodes = int(mesh_nodes_obj%size())
    ncells = int(mesh_cells_obj%size())
    ncells_shrinked = int(mesh_cells_obj%size())

    if (partit%mype == 0) then
      write(output_unit, '(A)') '  --> Setting up mesh from Atlas fesom-pi grid'
      write(output_unit, '(A,I0,A,I0)') '      Atlas mesh: ', ncells, ' cells, ', nnodes, ' nodes'
      write(output_unit, '(A)') '      Converting to FESOM format...'
    end if

    ! Convert atlas mesh to FESOM mesh
    call atlas_mesh_to_fesom_mesh(mesh_obj, zbar_default, partit, mesh)
    partit%myDim_elem2D_shrinked = ncells_shrinked
    deallocate(zbar_default)

    ! Store mesh for later use by compute_tracer_stats_atlas
    atlas_mesh_global = mesh_obj
    atlas_mesh_active = .true.
    call meshgen_obj%FINAL()
    call grid_obj%FINAL()
    ! Do NOT finalize mesh_obj - it's stored in atlas_mesh_global
#else
    call mesh_setup(partit, mesh)
#endif

  end subroutine mesh_setup_with_atlas

#ifdef ENABLE_ATLAS
  !> @brief Convert an Atlas mesh to a FESOM t_mesh
  !! @details Extracts node coordinates and triangle connectivity, asks Atlas
  !!          to build edge topology, maps it into the FESOM t_mesh arrays,
  !!          builds level and partition arrays, and calls the
  !!          standard mesh computation routines (mesh_areas,
  !!          mesh_auxiliary_arrays).
  !! @param[inout] atlas_msh  Atlas mesh object (RegularLonLat or similar)
  !! @param[in]    zbar       Vertical interface depths
  !! @param[inout] partit     FESOM partition structure (filled for npes=1)
  !! @param[inout] mesh3      Output FESOM t_mesh (must be uninitialised)
  subroutine atlas_mesh_to_fesom_mesh(atlas_msh, zbar, partit, mesh3)
    type(atlas_Mesh),   intent(inout) :: atlas_msh
    real(kind=MP),      intent(in)    :: zbar(:)
    type(t_partit),     intent(inout) :: partit
    type(t_mesh),       intent(inout) :: mesh3

    ! Atlas helper objects
    type(atlas_mesh_Nodes)             :: atlas_nodes
    type(atlas_mesh_Cells)             :: atlas_cells
    type(atlas_mesh_Edges)             :: atlas_edges
    type(atlas_MultiBlockConnectivity) :: atlas_edge_nodes
    type(atlas_MultiBlockConnectivity) :: atlas_edge_cells
    type(atlas_Field)                  :: node_global_index_field
    type(atlas_Field)                  :: cell_global_index_field
    type(atlas_Field)                  :: ghost_field
    integer(ATLAS_KIND_GIDX), pointer  :: node_global_index(:)
    integer(ATLAS_KIND_GIDX), pointer  :: cell_global_index(:)
    integer(c_int), pointer            :: node_ghost(:)
    integer(ATLAS_KIND_IDX), pointer   :: edge_nodes(:,:), edge_node_cols(:)
    integer(ATLAS_KIND_IDX), pointer   :: edge_cells(:,:), edge_cell_cols(:)

    ! Mesh dimensions
    integer :: nnodes, ncells, nl, edge2D_local, edge2D_in_local
    integer :: global_nnodes, global_ncells, local_global_nnodes
    integer :: local_global_ncells, file_unit, io_stat
    integer :: canonical_node_count
    integer, allocatable :: global_node_levels(:), global_cell_levels(:)
    integer, allocatable :: canonical_elements(:,:)
    real(kind=MP), allocatable :: canonical_lonlat(:,:)

    integer :: i, j, canonical_cell_count, global_node_id, canonical_node_flag
    integer :: tmp_el

    ! Edge classification
    integer :: n_internal
    integer, allocatable :: aux(:)
    integer :: e, n, k, q, e1
    integer :: elnodes(3), eledges(3)
    real(kind=WP) :: geographic_lon, geographic_lat, rotated_lon, rotated_lat

    ! Owned vs ghost node split
    integer :: n_owned, n_ghost, ierror

    real(MP), parameter :: deg2rad = acos(-1.0_MP) / 180.0_MP

    ! ------------------------------------------------------------------
    ! 1. Basic Atlas mesh dimensions
    ! ------------------------------------------------------------------
    atlas_nodes = atlas_msh%nodes()
    atlas_cells = atlas_msh%cells()
    nnodes = int(atlas_nodes%size())
    ncells = int(atlas_cells%size())
    nl = size(zbar)

    node_global_index_field = atlas_nodes%global_index()
    cell_global_index_field = atlas_cells%global_index()
    call node_global_index_field%data(node_global_index)
    call cell_global_index_field%data(cell_global_index)

    ghost_field = atlas_nodes%ghost()
    n_owned = 0
    n_ghost = 0
    call ghost_field%data(node_ghost)
    do n=1, nnodes
      if (node_ghost(n) == 0_c_int) then
        n_owned = n_owned + 1
        if (n_ghost > 0) then
          STOP 'Atlas mesh has owned nodes after ghost nodes; fesom_mesh requires all owned nodes first'
        end if
      else
        n_ghost = n_ghost + 1
      end if
    end do
    call ghost_field%final()
    if (partit%mype == 0) then
      write(output_unit, '(A,I0,A,I0)') &
        '  atlas_mesh_to_fesom_mesh: nnodes=', nnodes, ', ncells=', ncells
    end if

    ! ------------------------------------------------------------------
    ! 2. Canonical node coordinates: lon/lat degrees -> radians
    ! ------------------------------------------------------------------
    open(newunit=file_unit, file=trim(MeshPath)//'nod2d.out', &
         status='old', action='read', iostat=io_stat)
    if (io_stat /= 0) error stop 'Cannot open nod2d.out for Atlas mesh'
    read(file_unit, *, iostat=io_stat) canonical_node_count
    if (io_stat /= 0) error stop 'Cannot read nod2d.out header for Atlas mesh'
    allocate(canonical_lonlat(2, canonical_node_count))
    do n = 1, canonical_node_count
      read(file_unit, *, iostat=io_stat) i, canonical_lonlat(:, n), canonical_node_flag
      if (io_stat /= 0 .or. i /= n) error stop 'Cannot read nod2d.out for Atlas mesh'
    end do
    close(file_unit)

    mesh3%nod2D = nnodes
    allocate(mesh3%coord_nod2D(2, nnodes))
    call set_mesh_transform_matrix()
    do n = 1, nnodes
      geographic_lon = real(canonical_lonlat(1, int(node_global_index(n))) * deg2rad, WP)
      geographic_lat = real(canonical_lonlat(2, int(node_global_index(n))) * deg2rad, WP)
      if (force_rotation) then
        call g2r(geographic_lon, geographic_lat, rotated_lon, rotated_lat)
        mesh3%coord_nod2D(1, n) = real(rotated_lon, MP)
        mesh3%coord_nod2D(2, n) = real(rotated_lat, MP)
      else
        mesh3%coord_nod2D(1, n) = real(geographic_lon, MP)
        mesh3%coord_nod2D(2, n) = real(geographic_lat, MP)
      end if
    end do
    deallocate(canonical_lonlat)

    ! ------------------------------------------------------------------
    ! 3. Cell connectivity in canonical FESOM vertex order
    ! ------------------------------------------------------------------
    open(newunit=file_unit, file=trim(MeshPath)//'elem2d.out', &
         status='old', action='read', iostat=io_stat)
    if (io_stat /= 0) error stop 'Cannot open elem2d.out for Atlas mesh'
    read(file_unit, *, iostat=io_stat) canonical_cell_count
    if (io_stat /= 0) error stop 'Cannot read elem2d.out header for Atlas mesh'
    allocate(canonical_elements(3, canonical_cell_count))
    do e = 1, canonical_cell_count
      read(file_unit, *, iostat=io_stat) canonical_elements(:, e)
      if (io_stat /= 0) error stop 'Cannot read elem2d.out for Atlas mesh'
    end do
    close(file_unit)

    mesh3%elem2D = ncells
    allocate(mesh3%elem2D_nodes(3, ncells))
    do e = 1, ncells
      do k = 1, 3
        global_node_id = canonical_elements(k, int(cell_global_index(e)))
        mesh3%elem2D_nodes(k, e) = 0
        do n = 1, nnodes
          if (int(node_global_index(n)) == global_node_id) then
            mesh3%elem2D_nodes(k, e) = n
            exit
          end if
        end do
        if (mesh3%elem2D_nodes(k, e) == 0) then
          error stop 'Canonical FESOM cell references a node absent from Atlas mesh'
        end if
      end do
    end do
    deallocate(canonical_elements)

    ! ------------------------------------------------------------------
    ! 4. Vertical structure from aux3d.out, with a flat bottom
    ! ------------------------------------------------------------------
    mesh3%nl = nl
    allocate(mesh3%zbar(nl))
    mesh3%zbar = zbar
    allocate(mesh3%Z(nl-1))
    mesh3%Z = 0.5_MP * (mesh3%zbar(1:nl-1) + mesh3%zbar(2:nl))
    allocate(mesh3%depth(nnodes))
    mesh3%depth = mesh3%zbar(nl)
    allocate(mesh3%elem_depth(ncells))
    mesh3%elem_depth = mesh3%zbar(nl)

    ! ------------------------------------------------------------------
    ! 5. Build edge topology with Atlas
    ! ------------------------------------------------------------------
    call atlas_build_edges(atlas_msh)
    atlas_edges = atlas_msh%edges()
    atlas_edge_nodes = atlas_edges%node_connectivity()
    atlas_edge_cells = atlas_edges%cell_connectivity()
    call atlas_edge_nodes%padded_data(edge_nodes, edge_node_cols)
    call atlas_edge_cells%padded_data(edge_cells, edge_cell_cols)

    edge2D_local = int(atlas_edges%size())
    n_internal = 0
    do e = 1, edge2D_local
      if (edge_cell_cols(e) == 2 .and. edge_cells(2, e) > 0) then
        n_internal = n_internal + 1
      end if
    end do
    edge2D_in_local = n_internal

    ! Preserve Atlas-local edge ordering and orientation.
    allocate(mesh3%edges(2, edge2D_local))
    allocate(mesh3%edge_tri(2, edge2D_local))
    allocate(mesh3%elem_edges(3, ncells))
    mesh3%elem_edges = 0

    do e = 1, edge2D_local
      mesh3%edges(:, e) = int(edge_nodes(:, e))
      mesh3%edge_tri(1, e) = int(edge_cells(1, e))
      if (edge_cell_cols(e) == 2) then
        mesh3%edge_tri(2, e) = int(edge_cells(2, e))
      else
        mesh3%edge_tri(2, e) = 0
      end if
    end do

    mesh3%edge2D    = edge2D_local
    mesh3%edge2D_in = edge2D_in_local

    ! ------------------------------------------------------------------
    ! 6. Build elem_edges: elem_edges(q,e) = edge opposite to node q
    ! ------------------------------------------------------------------
    allocate(aux(ncells)); aux = 0
    do n = 1, edge2D_local
      do k = 1, 2
        q = mesh3%edge_tri(k, n)
        if (q > 0 .and. q <= ncells) then
          aux(q) = aux(q) + 1
          mesh3%elem_edges(aux(q), q) = n
        end if
      end do
    end do
    deallocate(aux)
    ! Reorder so elem_edges(q,e) does not contain node q
    do e = 1, ncells
      elnodes = mesh3%elem2D_nodes(:, e)
      eledges = mesh3%elem_edges(:, e)
      do q = 1, 3
        do k = 1, 3
          if (mesh3%edges(1,eledges(k)) /= elnodes(q) .and. &
              mesh3%edges(2,eledges(k)) /= elnodes(q)) then
            mesh3%elem_edges(q, e) = eledges(k); exit
          end if
        end do
      end do
    end do

    ! ------------------------------------------------------------------
    ! 7. Partition arrays: Atlas stores owned nodes before ghost nodes.
    ! ------------------------------------------------------------------
    partit%myDim_nod2D  = n_owned;      partit%eDim_nod2D   = n_ghost
    partit%myDim_elem2D = ncells;       partit%eDim_elem2D  = 0
    partit%eXDim_elem2D = 0
    partit%myDim_edge2D = edge2D_local; partit%eDim_edge2D  = 0

    allocate(partit%myList_nod2D(n_owned))
    do n = 1, n_owned
      partit%myList_nod2D(n) = int(node_global_index(n))
    end do
    allocate(partit%myList_elem2D(ncells))
    do n = 1, ncells
      partit%myList_elem2D(n) = int(cell_global_index(n))
    end do
    allocate(partit%myList_edge2D(edge2D_local))
    do n = 1, edge2D_local
      partit%myList_edge2D(n) = n
    end do
    ! CSR partition vector: size npes+1 with a balanced block distribution
    allocate(partit%part(partit%npes+1))
    partit%part(1) = 1
    do n = 1, partit%npes
      partit%part(n+1) = 1 + (n * nnodes) / partit%npes
    end do
    partit%part(partit%npes+1) = nnodes + 1

    ! Empty halo communication structures (npes=1: no halo exchange)
    partit%com_nod2D%rPEnum  = 0; partit%com_nod2D%sPEnum  = 0
    partit%com_elem2D%rPEnum = 0; partit%com_elem2D%sPEnum = 0
    partit%com_elem2D_full%rPEnum = 0; partit%com_elem2D_full%sPEnum = 0
    allocate(partit%com_nod2D%rlist(0))
    allocate(partit%com_nod2D%slist(0))
    allocate(partit%com_elem2D%rlist(0))
    allocate(partit%com_elem2D%slist(0))
    allocate(partit%com_elem2D_full%rlist(0))
    allocate(partit%com_elem2D_full%slist(0))
    allocate(partit%myInd_elem2D_shrinked(ncells))
    do n = 1, ncells; partit%myInd_elem2D_shrinked(n) = n; end do
    partit%myDim_elem2D_shrinked = ncells

    ! ------------------------------------------------------------------
    ! 8. Level arrays
    ! ------------------------------------------------------------------
    allocate(mesh3%ulevels(ncells));          mesh3%ulevels = 1
    allocate(mesh3%nlevels(ncells))
    allocate(mesh3%ulevels_nod2D(nnodes));    mesh3%ulevels_nod2D = 1
    allocate(mesh3%nlevels_nod2D(nnodes))
    allocate(mesh3%bc_index_nod2D(nnodes));   mesh3%bc_index_nod2D = 0

    local_global_nnodes = int(maxval(node_global_index))
    local_global_ncells = int(maxval(cell_global_index))
    call MPI_Allreduce(local_global_nnodes, global_nnodes, 1, MPI_INTEGER, &
                       MPI_MAX, partit%MPI_COMM_FESOM, ierror)
    call MPI_Allreduce(local_global_ncells, global_ncells, 1, MPI_INTEGER, &
                       MPI_MAX, partit%MPI_COMM_FESOM, ierror)
    allocate(global_node_levels(global_nnodes))
    allocate(global_cell_levels(global_ncells))

    if (partit%mype == 0) then
      open(newunit=file_unit, file=trim(MeshPath)//'nlvls.out', status='old', &
           action='read', iostat=io_stat)
      if (io_stat /= 0) error stop 'Cannot open nlvls.out for Atlas mesh'
      read(file_unit, *, iostat=io_stat) global_node_levels
      close(file_unit)
      if (io_stat /= 0) error stop 'Cannot read nlvls.out for Atlas mesh'

      open(newunit=file_unit, file=trim(MeshPath)//'elvls.out', status='old', &
           action='read', iostat=io_stat)
      if (io_stat /= 0) error stop 'Cannot open elvls.out for Atlas mesh'
      read(file_unit, *, iostat=io_stat) global_cell_levels
      close(file_unit)
      if (io_stat /= 0) error stop 'Cannot read elvls.out for Atlas mesh'
    end if
    call MPI_BCast(global_node_levels, global_nnodes, MPI_INTEGER, 0, &
                   partit%MPI_COMM_FESOM, ierror)
    call MPI_BCast(global_cell_levels, global_ncells, MPI_INTEGER, 0, &
                   partit%MPI_COMM_FESOM, ierror)

    do n = 1, nnodes
      mesh3%nlevels_nod2D(n) = global_node_levels(int(node_global_index(n)))
    end do
    do e = 1, ncells
      mesh3%nlevels(e) = global_cell_levels(int(cell_global_index(e)))
    end do
    deallocate(global_node_levels, global_cell_levels)

    ! ------------------------------------------------------------------
    ! 9. elem_neighbors and nod_in_elem2D
    ! ------------------------------------------------------------------
    allocate(mesh3%elem_neighbors(3, ncells))
    mesh3%elem_neighbors = 0
    do e = 1, ncells
      eledges = mesh3%elem_edges(:, e)
      do k = 1, 3
        e1 = mesh3%edge_tri(1, eledges(k))
        if (e1 == e) e1 = mesh3%edge_tri(2, eledges(k))
        mesh3%elem_neighbors(k, e) = e1
      end do
    end do

    allocate(mesh3%nod_in_elem2D_num(nnodes))
    mesh3%nod_in_elem2D_num = 0
    do e = 1, ncells
      do k = 1, 3
        n = mesh3%elem2D_nodes(k, e)
        mesh3%nod_in_elem2D_num(n) = mesh3%nod_in_elem2D_num(n) + 1
      end do
    end do
    allocate(mesh3%nod_in_elem2D(maxval(mesh3%nod_in_elem2D_num), nnodes))
    mesh3%nod_in_elem2D = 0
    mesh3%nod_in_elem2D_num = 0
    do e = 1, ncells
      do k = 1, 3
        n = mesh3%elem2D_nodes(k, e)
        mesh3%nod_in_elem2D_num(n) = mesh3%nod_in_elem2D_num(n) + 1
        mesh3%nod_in_elem2D(mesh3%nod_in_elem2D_num(n), n) = e
      end do
    end do

    do n = 1, nnodes
      do i = 2, mesh3%nod_in_elem2D_num(n)
        tmp_el = mesh3%nod_in_elem2D(i, n)
        j = i - 1
        do while (j >= 1 .and. &
                  cell_global_index(mesh3%nod_in_elem2D(j, n)) > cell_global_index(tmp_el))
          mesh3%nod_in_elem2D(j+1, n) = mesh3%nod_in_elem2D(j, n)
          j = j - 1
        end do
        mesh3%nod_in_elem2D(j+1, n) = tmp_el
      end do
    end do

    call node_global_index_field%final()
    call cell_global_index_field%final()

    ! ------------------------------------------------------------------
    ! 10. Mesh computation routines
    ! ------------------------------------------------------------------
    call init_mpi_types(partit, mesh3)
    call init_gatherLists(partit)
    call test_tri(partit, mesh3)
    call find_levels_min_e2n(partit, mesh3)
    call mesh_areas(partit, mesh3)
    call mesh_auxiliary_arrays(partit, mesh3)

    if (partit%mype == 0) then
      write(output_unit, '(A,I0,A,I0,A,I0)') &
        '  atlas_mesh_to_fesom_mesh complete: ', nnodes, ' nodes, ', &
        ncells, ' cells, ', edge2D_local, ' edges'
    end if

  end subroutine atlas_mesh_to_fesom_mesh
#endif

  !> @brief Compute tracer min/max/sum using Atlas NodeColumns
  !! @details Atlas mesh overlap nodes are marked as ghosts during conversion,
  !!          so NodeColumns reductions include each global node exactly once.
  !! @param[in] tracer_data 2D tracer array (nz, npoints)
  !! @param[in] n_owned Number of owned nodes on this rank
  !! @param[in] partit Partition structure (for MPI communicator)
  !! @param[out] tmin Minimum value (real(8))
  !! @param[out] tmax Maximum value (real(8))
  !! @param[out] tsum Sum of all values (real(8))
  subroutine compute_tracer_stats_atlas(tracer_data, n_owned, partit, tmin, tmax, tsum)
    real(kind=WP), intent(in) :: tracer_data(:,:)
    integer, intent(in) :: n_owned
    type(t_partit), intent(in) :: partit
    real(8), intent(out) :: tmin, tmax, tsum

#ifdef ENABLE_ATLAS
    type(atlas_functionspace_NodeColumns) :: fs
    type(atlas_Field) :: tracer_field
    real(WP), pointer :: atlas_tracer(:,:)
    real(WP) :: tmin_wp, tmax_wp, tsum_wp
    integer :: node_count
    integer :: ierr
    real(8) :: tmin_loc, tmax_loc, tsum_loc

    if (.not. atlas_fesom_active()) then
      tmin_loc = dble(minval(tracer_data(:, 1:n_owned)))
      tmax_loc = dble(maxval(tracer_data(:, 1:n_owned)))
      tsum_loc = sum(dble(tracer_data(:, 1:n_owned)))
      call MPI_Allreduce(tmin_loc, tmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                         partit%MPI_COMM_FESOM, ierr)
      call MPI_Allreduce(tmax_loc, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                         partit%MPI_COMM_FESOM, ierr)
      call MPI_Allreduce(tsum_loc, tsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                         partit%MPI_COMM_FESOM, ierr)
      return
    end if

    fs = atlas_functionspace_NodeColumns(atlas_mesh_global,halo=2)
    tracer_field = fs%create_field(name='tracer', kind=atlas_real(WP), &
                                   levels=size(tracer_data, 1))
    call tracer_field%data(atlas_tracer)
    node_count = min(size(atlas_tracer, 2), size(tracer_data, 2))
    if (n_owned > node_count) then
      error stop 'compute_tracer_stats_atlas: owned-node count exceeds field size'
    end if

    atlas_tracer = 0.0_WP
    atlas_tracer(:, 1:n_owned) = tracer_data(:, 1:n_owned)

    call fs%minimum(tracer_field, tmin_wp)
    call fs%maximum(tracer_field, tmax_wp)
    call fs%order_independent_sum(tracer_field, tsum_wp)
    tmin = dble(tmin_wp)
    tmax = dble(tmax_wp)
    tsum = dble(tsum_wp)

    call tracer_field%final()
    call fs%final()
#else
    integer :: ierr
    real(8) :: tmin_loc, tmax_loc, tsum_loc
    
    ! Compute local min/max/sum on this rank
    tmin_loc = dble(minval(tracer_data(:, 1:n_owned)))
    tmax_loc = dble(maxval(tracer_data(:, 1:n_owned)))
    tsum_loc = sum(dble(tracer_data(:, 1:n_owned)))
    
    ! Reduce across all ranks
    call MPI_Allreduce(tmin_loc, tmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                       partit%MPI_COMM_FESOM, ierr)
    call MPI_Allreduce(tmax_loc, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       partit%MPI_COMM_FESOM, ierr)
    call MPI_Allreduce(tsum_loc, tsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       partit%MPI_COMM_FESOM, ierr)
#endif
  end subroutine compute_tracer_stats_atlas

#ifdef ENABLE_ATLAS
  !> @brief Compute statistics directly on an Atlas field using NodeColumns
  !! @details Calls fs%minimum, fs%maximum, and an order-independent sum.
  !!          The field must be real(8) precision. Results are globally reduced
  !!          across all MPI ranks by Atlas internally.
  !! @param[inout] field  Atlas field to compute statistics on
  !! @param[in]    fs     NodeColumns function space owning the field
  !! @param[out]   tmin   Global minimum value
  !! @param[out]   tmax   Global maximum value
  !! @param[out]   tsum   Global sum of all owned-node values
  subroutine compute_field_stats_atlas(field, fs, tmin, tmax, tsum)
    type(atlas_Field),                     intent(inout) :: field
    type(atlas_functionspace_NodeColumns), intent(in)    :: fs
    real(8), intent(out) :: tmin, tmax, tsum
    call fs%minimum(field, tmin)
    call fs%maximum(field, tmax)
    call fs%order_independent_sum(field, tsum)
  end subroutine compute_field_stats_atlas
#endif

end module atlas_fesom_mesh_module
