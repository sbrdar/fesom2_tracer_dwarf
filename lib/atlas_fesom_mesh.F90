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
  use, intrinsic :: iso_c_binding, only: c_double, c_int
#endif

  implicit none
  private

  public :: mesh_setup_with_atlas, compute_tracer_stats_atlas
#ifdef ENABLE_ATLAS
  public :: atlas_mesh_to_fesom_mesh, compute_field_stats_atlas
  ! Module-level Atlas mesh (persists for stats computation)
  type(atlas_Mesh), save :: atlas_mesh_global
  logical, allocatable, save :: atlas_node_owned(:)
#endif

contains

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
    character(len=64)  :: env_value
    character(len=32)  :: use_fesom_dist_str
    character(len=10)  :: npes_string
    character(len=256) :: dist_mesh_dir, rpart_file
    integer :: ncells, nnodes
    integer :: nl_default
    real(kind=MP) :: max_depth_default
    type(atlas_mesh_Nodes) :: mesh_nodes_obj
    type(atlas_mesh_Cells) :: mesh_cells_obj
    ! Partition file variables
    integer :: funit, r, node_idx, npes_in_file, error_status, ierror
    integer, allocatable :: part_csr(:), part_per_node(:)
    logical :: have_dist
    integer :: n_owned_from_dist   ! owned-node count derived from part_csr for this rank

    ! Use sensible defaults for fesom-pi grid
    nl_default = 10
    max_depth_default = 1000.0_MP

    ! Read nl from aux3d.out, matching the non-atlas read_mesh path
    if (partit%mype == 0) then
      open(newunit=funit, file=trim(MeshPath)//'aux3d.out', status='old', &
           action='read', iostat=io_stat)
      if (io_stat == 0) then
        read(funit, *, iostat=io_stat) nl_default
        close(funit)
        if (io_stat /= 0 .or. nl_default < 3) nl_default = 10
      end if
    end if
    if (partit%npes > 1) &
      call MPI_BCast(nl_default, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)

    ! Initialize Atlas (loads plugins including atlas-fesom, registers grids)
    call atlas_initialize()

    ! Attempt to create fesom-pi grid; fall back to mesh_setup if unavailable
    if (partit%mype == 0) then
      write(output_unit, '(A)') '  --> Attempting to use Atlas fesom-pi grid...'
    end if

    grid_obj = atlas_Grid("fesom-pi")
    write(output_unit, *) '  --> Atlas fesom-pi grid found: ', grid_obj%name()

    ! Optionally drive Atlas mesh distribution from the standard rpart.out file.
    ! Set ATLAS_USE_FESOM_DIST=1 to enable; the file is read from the same
    ! dist_<npes>/ directory as the non-Atlas mesh_setup() path.
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
        ! Convert CSR (1-based offsets) to per-point 0-based partition indices
        ! part_csr(r)..part_csr(r+1)-1 are node global indices for rank r-1
        allocate(part_per_node(int(grid_obj%size())))
        do r = 1, partit%npes
          do node_idx = part_csr(r), part_csr(r+1) - 1
            part_per_node(node_idx) = r - 1   ! 0-based partition index for Atlas
          end do
        end do
        ! owned-node count for this rank: part_csr(mype+1) .. part_csr(mype+2)-1
        n_owned_from_dist = part_csr(partit%mype+2) - part_csr(partit%mype+1)
        deallocate(part_csr)
        dist_obj   = atlas_GridDistribution(part_per_node)
        deallocate(part_per_node)
        have_dist  = .true.
        if (partit%mype == 0) write(output_unit, '(3A,I0,A)') &
          '  --> ATLAS_USE_FESOM_DIST: using partition from ', &
          trim(rpart_file), ' (', partit%npes, ' partitions)'
      else
        if (allocated(part_csr)) deallocate(part_csr)
        n_owned_from_dist = 0
      end if
    else
      n_owned_from_dist = 0
    end if

    ! Create mesh generator and generate mesh
    meshgen_obj = atlas_MeshGenerator("fesom")
    if (have_dist) then
      mesh_obj = meshgen_obj%GENERATE(grid_obj, dist_obj)
      call dist_obj%final()
    else
      mesh_obj = meshgen_obj%GENERATE(grid_obj)
    end if

    ! Get mesh dimensions for logging
    mesh_nodes_obj = mesh_obj%nodes()
    mesh_cells_obj = mesh_obj%cells()
    nnodes = int(mesh_nodes_obj%size())
    ncells = int(mesh_cells_obj%size())

    if (partit%mype == 0) then
      write(output_unit, '(A)') '  --> Setting up mesh from Atlas fesom-pi grid'
      write(output_unit, '(A,I0,A,I0)') '      Atlas mesh: ', ncells, ' cells, ', nnodes, ' nodes'
      write(output_unit, '(A)') '      Converting to FESOM format...'
    end if

    ! Convert atlas mesh to FESOM mesh
    call atlas_mesh_to_fesom_mesh(mesh_obj, nl_default, max_depth_default, &
                                  n_owned_from_dist, partit, mesh)

    ! Store mesh for later use by compute_tracer_stats_atlas
    atlas_mesh_global = mesh_obj
    call meshgen_obj%FINAL()
    call grid_obj%FINAL()
    ! Do NOT finalize mesh_obj - it's stored in atlas_mesh_global
#else
    call mesh_setup(partit, mesh)
#endif

  end subroutine mesh_setup_with_atlas

#ifdef ENABLE_ATLAS
  !> @brief Convert an Atlas mesh to a FESOM t_mesh
  !! @details Extracts node coordinates (lon/lat) and triangle connectivity
  !!          from atlas_msh, builds the full FESOM t_mesh topology (edges,
  !!          elem_edges, level arrays, partition arrays) and calls the
  !!          standard mesh computation routines (mesh_areas,
  !!          mesh_auxiliary_arrays).
  !! @param[inout] atlas_msh  Atlas mesh object (RegularLonLat or similar)
  !! @param[in]    nl         Number of vertical levels (interfaces)
  !! @param[in]    max_depth  Flat-bottom depth (positive metres)
  !! @param[inout] partit     FESOM partition structure (filled for npes=1)
  !! @param[inout] mesh3      Output FESOM t_mesh (must be uninitialised)
  subroutine atlas_mesh_to_fesom_mesh(atlas_msh, nl, max_depth, n_owned_arg, partit, mesh3)
    type(atlas_Mesh),   intent(inout) :: atlas_msh
    integer,            intent(in)    :: nl
    real(kind=MP),      intent(in)    :: max_depth
    integer,            intent(in)    :: n_owned_arg  ! >0: owned count from dist; 0: single-rank
    type(t_partit),     intent(inout) :: partit
    type(t_mesh),       intent(inout) :: mesh3

    ! Atlas helper objects
    type(atlas_mesh_Nodes)             :: atlas_nodes
    type(atlas_mesh_Cells)             :: atlas_cells
    type(atlas_Field)                  :: lonlat_field
    type(atlas_Field)                  :: global_index_field
    type(atlas_MultiBlockConnectivity) :: atlas_conn
    real(c_double), pointer            :: lonlat(:,:)        ! (2, nnodes) degrees
    integer(ATLAS_KIND_GIDX), pointer  :: global_index(:)
    integer(ATLAS_KIND_IDX), pointer   :: conn_padded(:,:)   ! (3, ncells)  1-based
    integer(ATLAS_KIND_IDX), pointer   :: conn_ncols(:)      ! (ncells)

    ! Mesh dimensions
    integer :: nnodes, ncells, edge2D_local, edge2D_in_local

    ! Half-edge sorting arrays (3 per cell)
    integer :: n_hedges, i, j
    integer, allocatable :: he_lo(:), he_hi(:), he_el(:)
    integer :: tmp_lo, tmp_hi, tmp_el, tmp_lk, n_lo, n_hi

    ! Edge classification
    integer :: edge_count, n_internal, n_boundary, eidx
    logical, allocatable :: edge_is_boundary(:)
    integer, allocatable :: edge_perm(:), aux(:)
    integer :: e, n, k, q, e1
    integer :: elnodes(3), eledges(3)
    real(kind=WP) :: geographic_lon, geographic_lat, rotated_lon, rotated_lat

    ! Owned vs ghost node split
    integer :: n_owned, n_ghost, max_global_index, local_max_global_index, ierror
    integer, allocatable :: node_owner(:), global_node_owner(:)

    real(MP), parameter :: deg2rad = acos(-1.0_MP) / 180.0_MP

    ! ------------------------------------------------------------------
    ! 1. Basic Atlas mesh dimensions
    ! ------------------------------------------------------------------
    atlas_nodes = atlas_msh%nodes()
    atlas_cells = atlas_msh%cells()
    nnodes = int(atlas_nodes%size())
    ncells = int(atlas_cells%size())

    global_index_field = atlas_nodes%global_index()
    call global_index_field%data(global_index)
    local_max_global_index = int(maxval(global_index))
    call MPI_Allreduce(local_max_global_index, max_global_index, 1, MPI_INTEGER, &
                       MPI_MAX, partit%MPI_COMM_FESOM, ierror)
    allocate(node_owner(max_global_index), global_node_owner(max_global_index))
    node_owner = huge(0)
    do n = 1, nnodes
      node_owner(int(global_index(n))) = partit%mype
    end do
    call MPI_Allreduce(node_owner, global_node_owner, max_global_index, MPI_INTEGER, &
                       MPI_MIN, partit%MPI_COMM_FESOM, ierror)
    if (allocated(atlas_node_owned)) deallocate(atlas_node_owned)
    allocate(atlas_node_owned(nnodes))
    do n = 1, nnodes
      atlas_node_owned(n) = global_node_owner(int(global_index(n))) == partit%mype
    end do
    deallocate(node_owner, global_node_owner)
    call global_index_field%final()

    if (n_owned_arg > 0) then
      n_owned = n_owned_arg
    else
      n_owned = nnodes
    end if
    n_ghost = nnodes - n_owned

    if (partit%mype == 0) then
      write(output_unit, '(A,I0,A,I0)') &
        '  atlas_mesh_to_fesom_mesh: nnodes=', nnodes, ', ncells=', ncells
    end if

    ! ------------------------------------------------------------------
    ! 2. Node coordinates: lon/lat degrees -> radians -> coord_nod2D
    ! ------------------------------------------------------------------
    lonlat_field = atlas_nodes%lonlat()
    call lonlat_field%data(lonlat)   ! lonlat(1,n)=lon, lonlat(2,n)=lat [degrees]
    mesh3%nod2D = nnodes
    allocate(mesh3%coord_nod2D(2, nnodes))
    call set_mesh_transform_matrix()
    do n = 1, nnodes
      geographic_lon = real(real(lonlat(1, n), MP) * deg2rad, WP)
      geographic_lat = real(real(lonlat(2, n), MP) * deg2rad, WP)
      if (force_rotation) then
        call g2r(geographic_lon, geographic_lat, rotated_lon, rotated_lat)
        mesh3%coord_nod2D(1, n) = real(rotated_lon, MP)
        mesh3%coord_nod2D(2, n) = real(rotated_lat, MP)
      else
        mesh3%coord_nod2D(1, n) = real(geographic_lon, MP)
        mesh3%coord_nod2D(2, n) = real(geographic_lat, MP)
      end if
    end do
    call lonlat_field%final()

    ! ------------------------------------------------------------------
    ! 3. Cell connectivity: Atlas Fortran indices are already 1-based
    ! ------------------------------------------------------------------
    atlas_conn = atlas_cells%node_connectivity()
    call atlas_conn%padded_data(conn_padded, conn_ncols)
    ! conn_padded(k, e): k-th node of element e, 1-based; shape (3, ncells)
    mesh3%elem2D = ncells
    allocate(mesh3%elem2D_nodes(3, ncells))
    do e = 1, ncells
      do k = 1, 3
        mesh3%elem2D_nodes(k, e) = int(conn_padded(k, e))
      end do
    end do

    ! ------------------------------------------------------------------
    ! 4. Vertical structure (uniform, flat bottom)
    ! ------------------------------------------------------------------
    mesh3%nl = nl
    allocate(mesh3%zbar(nl))
    do k = 1, nl
      mesh3%zbar(k) = -real(k-1, MP) * max_depth / real(nl-1, MP)
    end do
    allocate(mesh3%Z(nl-1))
    do k = 1, nl-1
      mesh3%Z(k) = 0.5_MP * (mesh3%zbar(k) + mesh3%zbar(k+1))
    end do
    allocate(mesh3%depth(nnodes))
    mesh3%depth = -max_depth
    allocate(mesh3%elem_depth(ncells))
    mesh3%elem_depth = -max_depth

    ! ------------------------------------------------------------------
    ! 5. Build edges from element connectivity using half-edge sort
    !    Local edges of triangle (n1,n2,n3):
    !      k=1: (n1,n2)   k=2: (n2,n3)   k=3: (n1,n3)
    ! ------------------------------------------------------------------
    n_hedges = 3 * ncells
    allocate(he_lo(n_hedges), he_hi(n_hedges), he_el(n_hedges))

    do e = 1, ncells
      do k = 1, 3
        select case(k)
        case(1); n_lo = mesh3%elem2D_nodes(1,e); n_hi = mesh3%elem2D_nodes(2,e)
        case(2); n_lo = mesh3%elem2D_nodes(2,e); n_hi = mesh3%elem2D_nodes(3,e)
        case(3); n_lo = mesh3%elem2D_nodes(1,e); n_hi = mesh3%elem2D_nodes(3,e)
        end select
        i = (e-1)*3 + k
        he_lo(i) = min(n_lo, n_hi)
        he_hi(i) = max(n_lo, n_hi)
        he_el(i) = e
      end do
    end do

    ! Insertion sort on (he_lo, he_hi) — O(n^2), fine for test-case sizes
    do i = 2, n_hedges
      tmp_lo = he_lo(i); tmp_hi = he_hi(i); tmp_el = he_el(i)
      j = i - 1
      do while (j >= 1 .and. &
                (he_lo(j) > tmp_lo .or. &
                 (he_lo(j) == tmp_lo .and. he_hi(j) > tmp_hi)))
        he_lo(j+1) = he_lo(j); he_hi(j+1) = he_hi(j); he_el(j+1) = he_el(j)
        j = j - 1
      end do
      he_lo(j+1) = tmp_lo; he_hi(j+1) = tmp_hi; he_el(j+1) = tmp_el
    end do

    ! Count total / internal / boundary edges
    edge_count = 0; n_internal = 0; n_boundary = 0
    i = 1
    do while (i <= n_hedges)
      edge_count = edge_count + 1
      if (i < n_hedges .and. he_lo(i)==he_lo(i+1) .and. he_hi(i)==he_hi(i+1)) then
        n_internal = n_internal + 1; i = i + 2
      else
        n_boundary = n_boundary + 1; i = i + 1
      end if
    end do
    edge2D_local   = edge_count
    edge2D_in_local = n_internal

    ! Build permutation: internal edges first, boundary edges last
    allocate(edge_is_boundary(edge2D_local))
    edge_count = 0
    i = 1
    do while (i <= n_hedges)
      edge_count = edge_count + 1
      if (i < n_hedges .and. he_lo(i)==he_lo(i+1) .and. he_hi(i)==he_hi(i+1)) then
        edge_is_boundary(edge_count) = .false.; i = i + 2
      else
        edge_is_boundary(edge_count) = .true.;  i = i + 1
      end if
    end do
    allocate(edge_perm(edge2D_local))
    eidx = 0
    do n = 1, edge2D_local
      if (.not. edge_is_boundary(n)) then; eidx = eidx + 1; edge_perm(n) = eidx; end if
    end do
    do n = 1, edge2D_local
      if (edge_is_boundary(n))      then; eidx = eidx + 1; edge_perm(n) = eidx; end if
    end do

    ! Fill mesh3%edges and mesh3%edge_tri using permuted indices
    allocate(mesh3%edges(2, edge2D_local))
    allocate(mesh3%edge_tri(2, edge2D_local))
    allocate(mesh3%elem_edges(3, ncells))
    mesh3%elem_edges = 0

    edge_count = 0
    i = 1
    do while (i <= n_hedges)
      edge_count = edge_count + 1
      n = edge_perm(edge_count)
      if (i < n_hedges .and. he_lo(i)==he_lo(i+1) .and. he_hi(i)==he_hi(i+1)) then
        mesh3%edges(1,n) = he_lo(i);    mesh3%edges(2,n) = he_hi(i)
        mesh3%edge_tri(1,n) = he_el(i); mesh3%edge_tri(2,n) = he_el(i+1)
        i = i + 2
      else
        mesh3%edges(1,n) = he_lo(i);    mesh3%edges(2,n) = he_hi(i)
        mesh3%edge_tri(1,n) = he_el(i); mesh3%edge_tri(2,n) = 0
        i = i + 1
      end if
    end do
    deallocate(he_lo, he_hi, he_el, edge_is_boundary, edge_perm)

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
    do n = 1, n_owned; partit%myList_nod2D(n) = n; end do
    allocate(partit%myList_elem2D(ncells))
    do n = 1, ncells; partit%myList_elem2D(n) = n; end do
    allocate(partit%myList_edge2D(edge2D_local))
    do n = 1, edge2D_local; partit%myList_edge2D(n) = n; end do
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
    allocate(mesh3%nlevels(ncells));          mesh3%nlevels = nl
    allocate(mesh3%ulevels_nod2D(nnodes));    mesh3%ulevels_nod2D = 1
    allocate(mesh3%nlevels_nod2D(nnodes));    mesh3%nlevels_nod2D = nl
    allocate(mesh3%bc_index_nod2D(nnodes));   mesh3%bc_index_nod2D = 0

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

  !> @brief Compute tracer min/max/sum using MPI_Allreduce
  !! @details Atlas meshes can contain unmarked overlap nodes, so Atlas builds
  !!          reduce only the unique global-node ownership mask.
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
    integer :: ierr, node
    real(8) :: tmin_loc, tmax_loc, tsum_loc

    tmin_loc = huge(tmin_loc)
    tmax_loc = -huge(tmax_loc)
    tsum_loc = 0.0_8
    do node = 1, min(size(atlas_node_owned), size(tracer_data, 2))
      if (atlas_node_owned(node)) then
        tmin_loc = min(tmin_loc, dble(minval(tracer_data(:, node))))
        tmax_loc = max(tmax_loc, dble(maxval(tracer_data(:, node))))
        tsum_loc = tsum_loc + sum(dble(tracer_data(:, node)))
      end if
    end do
    call MPI_Allreduce(tmin_loc, tmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                       partit%MPI_COMM_FESOM, ierr)
    call MPI_Allreduce(tmax_loc, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       partit%MPI_COMM_FESOM, ierr)
    call MPI_Allreduce(tsum_loc, tsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       partit%MPI_COMM_FESOM, ierr)
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
  !! @details Calls fs%minimum, fs%maximum, fs%sum on the given field.
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
    call fs%sum(field, tsum)
  end subroutine compute_field_stats_atlas
#endif

end module atlas_fesom_mesh_module
