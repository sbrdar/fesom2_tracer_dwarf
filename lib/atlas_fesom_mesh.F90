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
  use o_PARAM, only: MP
  use oce_mesh_module, only: mesh_setup, test_tri, find_levels_min_e2n, &
                              mesh_areas, mesh_auxiliary_arrays
  use analytic_mesh_module, only: generate_analytic_mesh
  use g_rotate_grid, only: set_mesh_transform_matrix
  use par_support_interfaces, only: init_mpi_types, init_gatherLists
  use iso_fortran_env, only: output_unit

#ifdef ENABLE_ATLAS
  use atlas_module
  use, intrinsic :: iso_c_binding, only: c_double
#endif

  implicit none
  private

  public :: mesh_setup_with_atlas
#ifdef ENABLE_ATLAS
  public :: atlas_mesh_to_fesom_mesh
#endif

contains

  !> @brief Mesh setup entry point with optional Atlas path
  !! @param[inout] partit Partition structure
  !! @param[inout] mesh Mesh structure
  subroutine mesh_setup_with_atlas(partit, mesh)
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(inout), target :: mesh

#ifdef ENABLE_ATLAS
    type(atlas_StructuredGrid)  :: atlas_grid
    type(atlas_MeshGenerator)   :: atlas_meshgen
    type(atlas_Mesh)            :: atlas_mesh

    integer :: nx, ny, nl
    integer :: io_stat
    character(len=64)  :: env_value
    real(kind=MP), parameter :: Lx = 100000.0_MP
    real(kind=MP), parameter :: Ly = 100000.0_MP
    real(kind=MP), parameter :: max_depth = 1000.0_MP

    ! Defaults can be overridden via environment for quick experiments.
    nx = 20
    ny = 20
    nl = 10

    call get_environment_variable('FESOM_ATLAS_NX', env_value, status=io_stat)
    if (io_stat == 0) then
      read(env_value, *, iostat=io_stat) nx
    end if
    call get_environment_variable('FESOM_ATLAS_NY', env_value, status=io_stat)
    if (io_stat == 0) then
      read(env_value, *, iostat=io_stat) ny
    end if
    call get_environment_variable('FESOM_ATLAS_NL', env_value, status=io_stat)
    if (io_stat == 0) then
      read(env_value, *, iostat=io_stat) nl
    end if

    ! if (nx < 2 .or. ny < 2 .or. nl < 3) then
      if (partit%mype == 0) then
        write(output_unit, '(A)') 'Atlas mesh setup: invalid dimensions, falling back to mesh_setup()'
      end if
      call mesh_setup(partit, mesh)
      return
    ! end if
    if (partit%mype == 0) then
        write(output_unit, '(A,I0)') 'partit%part :', partit%part
    end if

    ! Create regular lon-lat grid
    atlas_grid = ATLAS_REGULARLONLATGRID(nx, ny)
    
    ! Create mesh generator and generate mesh
    atlas_meshgen = ATLAS_MESHGENERATOR("structured")
    atlas_mesh = atlas_meshgen%GENERATE(atlas_grid)

    if (partit%mype == 0) then
      write(output_unit, '(A)') '  --> Setting up mesh from Atlas generated grid'
      write(output_unit, '(A,I0,A,I0)') '      Atlas grid: ', nx, ' x ', ny
      write(output_unit, '(A)') '      Atlas mesh generated successfully'
    end if

    call generate_analytic_mesh(nx, ny, nl, Lx, Ly, max_depth, partit, mesh)

    ! call atlas_mesh%FINAL()
    call atlas_meshgen%FINAL()
    call atlas_grid%FINAL()
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
  subroutine atlas_mesh_to_fesom_mesh(atlas_msh, nl, max_depth, partit, mesh3)
    type(atlas_Mesh),   intent(inout) :: atlas_msh
    integer,            intent(in)    :: nl
    real(kind=MP),      intent(in)    :: max_depth
    type(t_partit),     intent(inout) :: partit
    type(t_mesh),       intent(inout) :: mesh3

    ! Atlas helper objects
    type(atlas_mesh_Nodes)             :: atlas_nodes
    type(atlas_mesh_Cells)             :: atlas_cells
    type(atlas_Field)                  :: atlas_field
    type(atlas_MultiBlockConnectivity) :: atlas_conn
    real(c_double), pointer            :: lonlat(:,:)        ! (2, nnodes) degrees
    integer(ATLAS_KIND_IDX), pointer   :: conn_padded(:,:)   ! (3, ncells)  0-based
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

    real(MP), parameter :: deg2rad = acos(-1.0_MP) / 180.0_MP

    ! ------------------------------------------------------------------
    ! 1. Basic Atlas mesh dimensions
    ! ------------------------------------------------------------------
    atlas_nodes = atlas_msh%nodes()
    atlas_cells = atlas_msh%cells()
    nnodes = int(atlas_nodes%size())
    ncells = int(atlas_cells%size())

    if (partit%mype == 0) then
      write(output_unit, '(A,I0,A,I0)') &
        '  atlas_mesh_to_fesom_mesh: nnodes=', nnodes, ', ncells=', ncells
    end if

    ! ------------------------------------------------------------------
    ! 2. Node coordinates: lon/lat degrees -> radians -> coord_nod2D
    ! ------------------------------------------------------------------
    atlas_field = atlas_nodes%lonlat()
    call atlas_field%data(lonlat)   ! lonlat(1,n)=lon, lonlat(2,n)=lat [degrees]
    mesh3%nod2D = nnodes
    allocate(mesh3%coord_nod2D(2, nnodes))
    do n = 1, nnodes
      mesh3%coord_nod2D(1, n) = real(lonlat(1, n), MP) * deg2rad
      mesh3%coord_nod2D(2, n) = real(lonlat(2, n), MP) * deg2rad
    end do
    call atlas_field%final()

    ! ------------------------------------------------------------------
    ! 3. Cell connectivity: 0-based Atlas indices -> 1-based elem2D_nodes
    ! ------------------------------------------------------------------
    atlas_conn = atlas_cells%node_connectivity()
    call atlas_conn%padded_data(conn_padded, conn_ncols)
    ! conn_padded(k, e): k-th node of element e, 0-based; shape (3, ncells)
    mesh3%elem2D = ncells
    allocate(mesh3%elem2D_nodes(3, ncells))
    do e = 1, ncells
      do k = 1, 3
        mesh3%elem2D_nodes(k, e) = int(conn_padded(k, e)) + 1
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
    ! 7. Partition arrays (npes=1: single rank owns all)
    ! ------------------------------------------------------------------
    partit%myDim_nod2D  = nnodes;       partit%eDim_nod2D   = 0
    partit%myDim_elem2D = ncells;       partit%eDim_elem2D  = 0
    partit%eXDim_elem2D = 0
    partit%myDim_edge2D = edge2D_local; partit%eDim_edge2D  = 0

    allocate(partit%myList_nod2D(nnodes))
    do n = 1, nnodes; partit%myList_nod2D(n) = n; end do
    allocate(partit%myList_elem2D(ncells))
    do n = 1, ncells; partit%myList_elem2D(n) = n; end do
    allocate(partit%myList_edge2D(edge2D_local))
    do n = 1, edge2D_local; partit%myList_edge2D(n) = n; end do
    allocate(partit%part(2))
    partit%part(1) = 1; partit%part(2) = nnodes + 1

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
    call set_mesh_transform_matrix()
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

end module atlas_fesom_mesh_module
