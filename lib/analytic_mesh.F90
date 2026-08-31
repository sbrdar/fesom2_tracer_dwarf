!> @file analytic_mesh.F90
!! @brief Analytic mesh generator for FESOM2 tracer dwarf
!! @details Generates a parametrized rectangular triangular mesh entirely
!!          in-memory, with no file I/O. Enables npes=1 testing with
!!          arbitrary mesh sizes (nx, ny, nl).

module analytic_mesh_module
  use MOD_MESH
  use MOD_PARTIT
  use MOD_PARSUP
  use oce_mesh_module
  use g_config
  use g_rotate_grid
  use g_comm_auto
  use par_support_interfaces
  use o_PARAM
  use mpi
#ifdef USE_HALF_PRECISION
  use hp_math_intrinsics
#endif

  implicit none
  private
  public :: generate_analytic_mesh

contains

  !> @brief Generate an analytic rectangular triangular mesh
  !! @details Builds a Cartesian rectangular domain triangulated into right
  !!          triangles, with uniform vertical levels and flat bottom.
  !!          Calls existing mesh computation routines for derived quantities.
  !!
  !! Mesh topology for grid points (i,j) with i=1..nx, j=1..ny:
  !!
  !!   (i,j+1)-----(i+1,j+1)     Node index: node(i,j) = (j-1)*nx + i
  !!     | \  upper  |
  !!     |   \       |            nod2D  = nx * ny
  !!     |     \     |            elem2D = 2 * (nx-1) * (ny-1)
  !!     |  lower \  |
  !!     |         \ |
  !!   (i,j)------(i+1,j)
  !!
  !! Lower triangle: nodes (i,j), (i+1,j), (i+1,j+1) -- CCW
  !! Upper triangle: nodes (i,j), (i+1,j+1), (i,j+1) -- CCW
  !!
  !! Boundary condition treatment (non-periodic):
  !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !! The rectangular domain has closed-wall (no-normal-flux) boundaries on all
  !! four sides. This is enforced implicitly through the mesh topology, not via
  !! explicit BC arrays. The mechanism has three parts:
  !!
  !! 1. edge_tri(2,n) == 0 for boundary edges:
  !!    Edges on the domain perimeter have only one adjacent element. The
  !!    advection driver (oce_adv_tra_driver.F90) checks `if(el(2)>0)` before
  !!    accessing the second element's level range. When el(2)==0, the flux
  !!    contribution from the missing side is zero.
  !!
  !! 2. edge_cross_dxdy(3:4,n) = 0 for boundary edges:
  !!    mesh_auxiliary_arrays (oce_mesh.F90) sets the geometric cross-distance
  !!    vector to zero for the missing element side. This zeros out gradient
  !!    contributions from the non-existent neighbor.
  !!
  !! 3. nboundary_lay array (oce_muscl_adv.F90):
  !!    During FCT/MUSCL initialization, nodes on boundary edges get
  !!    nboundary_lay=0, which disables higher-order reconstruction at those
  !!    nodes and falls back to first-order upwind (zero-gradient condition).
  !!
  !! This works for any nx, ny: boundary edges are wherever edge_tri has a
  !! zero entry, which is automatically the perimeter of the rectangle.
  !! Total boundary edges = 2*(nx-1) + 2*(ny-1).
  !!
  !! Periodic boundary condition treatment (periodic=.true.):
  !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !! Uses temporary halo (ghost) nodes with shifted coordinates during mesh
  !! construction, then remaps all connectivity to use real nodes after
  !! geometry is computed. This eliminates boundary artifacts.
  !!
  !! Real nodes: (nx-1)*(ny-1) — grid points i=1..nx-1, j=1..ny-1
  !! Right halos: ny-1 — copies of i=1 column with x += Lx
  !! Top halos: nx-1 — copies of j=1 row with y += Ly
  !! Corner halo: 1 — copy of (1,1) with both shifts
  !! All edges are internal (no boundary edges).
  !!
  !! After mesh_areas + mesh_auxiliary_arrays:
  !! - Fix edge_cross_dxdy with minimum-image convention
  !! - Merge halo node areas into real nodes
  !! - Remap elem2D_nodes and edges to real node indices
  !! - Set eDim_nod2D = 0
  !!
  !! IMPORTANT convention — edge_tri(1,n) must always be valid (>0):
  !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !! FESOM2's mesh_auxiliary_arrays (oce_mesh.F90, line ~2422) accesses
  !! center_x(edge_tri(1,n)) unconditionally, without a guard. Only
  !! edge_tri(2,n) is checked with `if(el(2)>0)`. This means boundary
  !! edges MUST have the valid element in position 1 and the zero in
  !! position 2. If this convention is violated (edge_tri(1,n)==0), you
  !! get an out-of-bounds access to center_x(0) which may silently
  !! corrupt memory on small meshes and segfault on larger ones.
  !!
  !! The edge generation below naturally produces some boundary edges
  !! with the zero in position 1 (e.g. bottom-row horizontal edges where
  !! e_upper=0). A swap step after edge generation enforces the convention
  !! by moving any zero from position 1 to position 2.
  !!
  !! @param[in] nx Number of grid points in x direction
  !! @param[in] ny Number of grid points in y direction
  !! @param[in] nl_in Number of vertical levels (must be >= 3)
  !! @param[in] Lx Domain size in x direction [meters]
  !! @param[in] Ly Domain size in y direction [meters]
  !! @param[in] max_depth Maximum depth [meters, positive]
  !! @param[inout] partit Partition structure (filled for npes=1)
  !! @param[inout] mesh Mesh structure (fully populated)
  !! @param[in] periodic Optional: if .true., use doubly-periodic boundaries
  subroutine generate_analytic_mesh(nx, ny, nl_in, Lx, Ly, max_depth, partit, mesh, periodic)
    use iso_fortran_env, only: output_unit

    integer,       intent(in)    :: nx, ny, nl_in
    real(kind=MP), intent(in)    :: Lx, Ly, max_depth
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(inout), target :: mesh
    logical,        intent(in), optional :: periodic

    integer :: i, j, k, n, e, eidx
    integer :: nod2D_local, elem2D_local, edge2D_local, edge2D_in_local
    real(kind=MP) :: dx, dy
    integer :: n_horiz, n_vert, n_diag, n_total_edges
    integer :: edge_count, n_internal, n_boundary
    integer, allocatable :: edge_nodes(:,:), edge_elems(:,:)
    logical, allocatable :: edge_is_boundary(:)
    integer, allocatable :: edge_perm(:)
    integer, allocatable :: edge_nodes_sorted(:,:), edge_elems_sorted(:,:)
    integer :: n1, n2, n3, n4, e_lower, e_upper
    integer, allocatable :: aux(:)
    integer :: elnodes(3), eledges(3), q
    integer :: elem, elem1, j_nb, node, max_nod_in_elem

    ! Periodic mesh variables
    logical :: is_periodic
    integer :: eDim, nod2D_total

    is_periodic = .false.
    if (present(periodic)) is_periodic = periodic

    if (partit%mype == 0) then
      write(output_unit, '(A)') ''
      write(output_unit, '(A)') '================================================'
      write(output_unit, '(A)') 'Generating analytic mesh'
      write(output_unit, '(A)') '================================================'
      write(output_unit, '(A,I6,A,I6,A,I4)') '  Grid: nx=', nx, ', ny=', ny, ', nl=', nl_in
      write(output_unit, '(A,E12.4,A,E12.4)') '  Domain: Lx=', Lx, ', Ly=', Ly
      write(output_unit, '(A,E12.4)') '  Max depth: ', max_depth
      if (is_periodic) then
        write(output_unit, '(A)') '  Boundary: DOUBLY-PERIODIC'
      else
        write(output_unit, '(A)') '  Boundary: CLOSED-WALL'
      end if
    end if

    ! ========================================
    ! 1a. Set global config flags
    ! ========================================
    cartesian = .true.
    force_rotation = .false.
    fplane = .true.
    use_cavity = .false.
    use_ice = .false.
    use_sw_pene = .false.
    use_depthfile = 'aux3d'
    toy_ocean = .true.
    flag_debug = .false.

    ! ========================================
    ! 1b. Compute mesh sizes
    ! ========================================
    if (is_periodic) then
      nod2D_local  = (nx-1) * (ny-1)
      eDim         = (nx-1) + (ny-1) + 1  ! right + top + corner halos
      n_horiz      = (nx-1) * (ny-1)
      n_vert       = (nx-1) * (ny-1)
      n_diag       = (nx-1) * (ny-1)
    else
      nod2D_local  = nx * ny
      eDim         = 0
      n_horiz      = (nx-1) * ny
      n_vert       = nx * (ny-1)
      n_diag       = (nx-1) * (ny-1)
    end if
    nod2D_total  = nod2D_local + eDim
    elem2D_local = 2 * (nx-1) * (ny-1)
    n_total_edges = n_horiz + n_vert + n_diag

    ! ========================================
    ! 1b. Fill t_partit for npes=1
    ! ========================================
    mesh%nod2D  = nod2D_local
    mesh%elem2D = elem2D_local

    partit%myDim_nod2D  = nod2D_local
    partit%eDim_nod2D   = eDim

    partit%myDim_elem2D = elem2D_local
    partit%eDim_elem2D  = 0
    partit%eXDim_elem2D = 0

    ! Edge dimensions (set after edge construction)
    partit%myDim_edge2D = n_total_edges
    partit%eDim_edge2D  = 0

    ! Identity mapping lists
    allocate(partit%myList_nod2D(nod2D_total))
    do n = 1, nod2D_total
      partit%myList_nod2D(n) = n
    end do

    allocate(partit%myList_elem2D(elem2D_local))
    do n = 1, elem2D_local
      partit%myList_elem2D(n) = n
    end do

    allocate(partit%myList_edge2D(n_total_edges))
    do n = 1, n_total_edges
      partit%myList_edge2D(n) = n
    end do

    ! Partition vector for npes=1
    allocate(partit%part(2))
    partit%part(1) = 1
    partit%part(2) = nod2D_local + 1

    ! Communication structures: no halo exchange needed for npes=1
    partit%com_nod2D%rPEnum = 0
    partit%com_nod2D%sPEnum = 0
    partit%com_elem2D%rPEnum = 0
    partit%com_elem2D%sPEnum = 0
    partit%com_elem2D_full%rPEnum = 0
    partit%com_elem2D_full%sPEnum = 0

    ! Allocate empty rlist/slist (some routines may reference these)
    allocate(partit%com_nod2D%rlist(0))
    allocate(partit%com_nod2D%slist(0))
    allocate(partit%com_elem2D%rlist(0))
    allocate(partit%com_elem2D%slist(0))
    allocate(partit%com_elem2D_full%rlist(0))
    allocate(partit%com_elem2D_full%slist(0))

    ! Shrunk elem list for npes=1: all elements
    partit%myDim_elem2D_shrinked = elem2D_local
    allocate(partit%myInd_elem2D_shrinked(elem2D_local))
    do n = 1, elem2D_local
      partit%myInd_elem2D_shrinked(n) = n
    end do

    ! Representative checksum
    mesh%representative_checksum = ''

    if (partit%mype == 0) then
      write(output_unit, '(A,I8)') '  nod2D  = ', nod2D_local
      if (is_periodic) then
        write(output_unit, '(A,I8)') '  halo nodes = ', eDim
      end if
      write(output_unit, '(A,I8)') '  elem2D = ', elem2D_local
      write(output_unit, '(A,I8)') '  edge2D = ', n_total_edges
    end if

    ! ========================================
    ! 1c. Generate node coordinates
    ! ========================================
    dx = Lx / real(nx - 1, MP)
    dy = Ly / real(ny - 1, MP)

    allocate(mesh%coord_nod2D(2, nod2D_total))

    if (is_periodic) then
      ! Real nodes: i=1..nx-1, j=1..ny-1
      do j = 1, ny-1
        do i = 1, nx-1
          n = (j-1)*(nx-1) + i
          mesh%coord_nod2D(1, n) = real(i-1, MP) * dx / r_earth
          mesh%coord_nod2D(2, n) = real(j-1, MP) * dy / r_earth
        end do
      end do
      ! Right halos: copies of i=1 column with x += Lx
      do j = 1, ny-1
        n = nod2D_local + j
        mesh%coord_nod2D(1, n) = Lx / r_earth
        mesh%coord_nod2D(2, n) = real(j-1, MP) * dy / r_earth
      end do
      ! Top halos: copies of j=1 row with y += Ly
      do i = 1, nx-1
        n = nod2D_local + (ny-1) + i
        mesh%coord_nod2D(1, n) = real(i-1, MP) * dx / r_earth
        mesh%coord_nod2D(2, n) = Ly / r_earth
      end do
      ! Corner halo: copy of (1,1) with both shifts
      mesh%coord_nod2D(1, nod2D_total) = Lx / r_earth
      mesh%coord_nod2D(2, nod2D_total) = Ly / r_earth
    else
      do j = 1, ny
        do i = 1, nx
          n = (j-1)*nx + i
          mesh%coord_nod2D(1, n) = real(i-1, MP) * dx / r_earth
          mesh%coord_nod2D(2, n) = real(j-1, MP) * dy / r_earth
        end do
      end do
    end if

    ! ========================================
    ! 1d. Generate element connectivity
    ! ========================================
    allocate(mesh%elem2D_nodes(3, elem2D_local))
    e = 0
    if (is_periodic) then
      do j = 1, ny-1
        do i = 1, nx-1
          n1 = pnode(i,   j)
          n2 = pnode(i+1, j)
          n3 = pnode(i+1, j+1)
          n4 = pnode(i,   j+1)

          ! Lower triangle (CW ordering): n1, n3, n2
          e = e + 1
          mesh%elem2D_nodes(1, e) = n1
          mesh%elem2D_nodes(2, e) = n3
          mesh%elem2D_nodes(3, e) = n2

          ! Upper triangle (CW ordering): n1, n4, n3
          e = e + 1
          mesh%elem2D_nodes(1, e) = n1
          mesh%elem2D_nodes(2, e) = n4
          mesh%elem2D_nodes(3, e) = n3
        end do
      end do
    else
      do j = 1, ny-1
        do i = 1, nx-1
          n1 = (j-1)*nx + i       ! (i,j)
          n2 = (j-1)*nx + i + 1   ! (i+1,j)
          n3 = j*nx + i + 1       ! (i+1,j+1)
          n4 = j*nx + i           ! (i,j+1)

          ! Lower triangle (CW ordering): n1, n3, n2
          e = e + 1
          mesh%elem2D_nodes(1, e) = n1
          mesh%elem2D_nodes(2, e) = n3
          mesh%elem2D_nodes(3, e) = n2

          ! Upper triangle (CW ordering): n1, n4, n3
          e = e + 1
          mesh%elem2D_nodes(1, e) = n1
          mesh%elem2D_nodes(2, e) = n4
          mesh%elem2D_nodes(3, e) = n3
        end do
      end do
    end if

    ! ========================================
    ! 1e. Generate vertical levels
    ! ========================================
    mesh%nl = nl_in

    allocate(mesh%zbar(nl_in))
    do k = 1, nl_in
      mesh%zbar(k) = -real(k-1, MP) * max_depth / real(nl_in - 1, MP)
    end do

    allocate(mesh%Z(nl_in - 1))
    do k = 1, nl_in - 1
      mesh%Z(k) = 0.5_MP * (mesh%zbar(k) + mesh%zbar(k+1))
    end do

    ! Depth at nodes (flat bottom, negative)
    allocate(mesh%depth(nod2D_total))
    mesh%depth(:) = -max_depth

    ! ========================================
    ! 1f. Generate edges and edge_tri
    ! ========================================
    allocate(edge_nodes(2, n_total_edges))
    allocate(edge_elems(2, n_total_edges))
    edge_elems = 0
    edge_count = 0

    if (is_periodic) then
      ! --- Periodic edge generation: all edges have 2 adjacent elements ---

      ! Horizontal edges: (i,j)-(i+1,j) for i=1..nx-1, j=1..ny-1
      do j = 1, ny-1
        do i = 1, nx-1
          edge_count = edge_count + 1
          edge_nodes(1, edge_count) = pnode(i, j)
          edge_nodes(2, edge_count) = pnode(i+1, j)

          ! Below: lower tri of quad(i,j)
          e_lower = 2*((j-1)*(nx-1) + (i-1)) + 1
          ! Above: upper tri of quad(i,j-1), wrapping j-1=0 to ny-1
          if (j >= 2) then
            e_upper = 2*((j-2)*(nx-1) + (i-1)) + 2
          else
            e_upper = 2*((ny-2)*(nx-1) + (i-1)) + 2
          end if

          edge_elems(1, edge_count) = e_upper
          edge_elems(2, edge_count) = e_lower
        end do
      end do

      ! Vertical edges: (i,j)-(i,j+1) for i=1..nx-1, j=1..ny-1
      do j = 1, ny-1
        do i = 1, nx-1
          edge_count = edge_count + 1
          edge_nodes(1, edge_count) = pnode(i, j)
          edge_nodes(2, edge_count) = pnode(i, j+1)

          ! Right: upper tri of quad(i,j)
          e_upper = 2*((j-1)*(nx-1) + (i-1)) + 2
          ! Left: lower tri of quad(i-1,j), wrapping i-1=0 to nx-1
          if (i >= 2) then
            e_lower = 2*((j-1)*(nx-1) + (i-2)) + 1
          else
            e_lower = 2*((j-1)*(nx-1) + (nx-2)) + 1
          end if

          edge_elems(1, edge_count) = e_upper
          edge_elems(2, edge_count) = e_lower
        end do
      end do

      ! Diagonal edges: (i,j)-(i+1,j+1) for i=1..nx-1, j=1..ny-1
      do j = 1, ny-1
        do i = 1, nx-1
          edge_count = edge_count + 1
          edge_nodes(1, edge_count) = pnode(i, j)
          edge_nodes(2, edge_count) = pnode(i+1, j+1)

          e_lower = 2*((j-1)*(nx-1) + (i-1)) + 1
          e_upper = 2*((j-1)*(nx-1) + (i-1)) + 2

          edge_elems(1, edge_count) = e_upper
          edge_elems(2, edge_count) = e_lower
        end do
      end do

      ! All edges are internal — no boundary sort needed
      edge2D_local = edge_count
      edge2D_in_local = edge_count
      n_internal = edge_count
      n_boundary = 0

      mesh%edge2D = edge2D_local
      mesh%edge2D_in = edge2D_in_local

      allocate(mesh%edges(2, edge2D_local))
      allocate(mesh%edge_tri(2, edge2D_local))
      do n = 1, edge2D_local
        mesh%edges(1, n) = edge_nodes(1, n)
        mesh%edges(2, n) = edge_nodes(2, n)
        mesh%edge_tri(1, n) = edge_elems(1, n)
        mesh%edge_tri(2, n) = edge_elems(2, n)
      end do

      deallocate(edge_nodes, edge_elems)

    else
      ! --- Non-periodic edge generation (original code) ---

      ! Horizontal edges: connect (i,j) -- (i+1,j) for i=1..nx-1, j=1..ny
      do j = 1, ny
        do i = 1, nx-1
          edge_count = edge_count + 1
          n1 = (j-1)*nx + i
          n2 = (j-1)*nx + i + 1
          edge_nodes(1, edge_count) = n1
          edge_nodes(2, edge_count) = n2

          if (j <= ny-1) then
            e_lower = 2*((j-1)*(nx-1) + (i-1)) + 1
          else
            e_lower = 0
          end if
          if (j >= 2) then
            e_upper = 2*((j-2)*(nx-1) + (i-1)) + 2
          else
            e_upper = 0
          end if

          edge_elems(1, edge_count) = e_upper
          edge_elems(2, edge_count) = e_lower
        end do
      end do

      ! Vertical edges: connect (i,j) -- (i,j+1) for i=1..nx, j=1..ny-1
      do j = 1, ny-1
        do i = 1, nx
          edge_count = edge_count + 1
          n1 = (j-1)*nx + i
          n2 = j*nx + i
          edge_nodes(1, edge_count) = n1
          edge_nodes(2, edge_count) = n2

          if (i <= nx-1) then
            e_upper = 2*((j-1)*(nx-1) + (i-1)) + 2
          else
            e_upper = 0
          end if
          if (i >= 2) then
            e_lower = 2*((j-1)*(nx-1) + (i-2)) + 1
          else
            e_lower = 0
          end if

          edge_elems(1, edge_count) = e_upper
          edge_elems(2, edge_count) = e_lower
        end do
      end do

      ! Diagonal edges: connect (i,j) -- (i+1,j+1) for i=1..nx-1, j=1..ny-1
      do j = 1, ny-1
        do i = 1, nx-1
          edge_count = edge_count + 1
          n1 = (j-1)*nx + i
          n2 = j*nx + i + 1
          edge_nodes(1, edge_count) = n1
          edge_nodes(2, edge_count) = n2

          e_lower = 2*((j-1)*(nx-1) + (i-1)) + 1
          e_upper = 2*((j-1)*(nx-1) + (i-1)) + 2

          edge_elems(1, edge_count) = e_upper
          edge_elems(2, edge_count) = e_lower
        end do
      end do

      ! Enforce FESOM2 convention: edge_tri(1,n) must always be valid (>0)
      do n = 1, n_total_edges
        if (edge_elems(1, n) == 0 .and. edge_elems(2, n) > 0) then
          edge_elems(1, n) = edge_elems(2, n)
          edge_elems(2, n) = 0
        end if
      end do

      ! Classify and sort: internal edges first, then boundary
      allocate(edge_is_boundary(n_total_edges))
      n_internal = 0
      n_boundary = 0
      do n = 1, n_total_edges
        if (edge_elems(1, n) > 0 .and. edge_elems(2, n) > 0) then
          edge_is_boundary(n) = .false.
          n_internal = n_internal + 1
        else
          edge_is_boundary(n) = .true.
          n_boundary = n_boundary + 1
        end if
      end do

      allocate(edge_perm(n_total_edges))
      edge2D_in_local = n_internal
      eidx = 0
      do n = 1, n_total_edges
        if (.not. edge_is_boundary(n)) then
          eidx = eidx + 1
          edge_perm(n) = eidx
        end if
      end do
      do n = 1, n_total_edges
        if (edge_is_boundary(n)) then
          eidx = eidx + 1
          edge_perm(n) = eidx
        end if
      end do

      edge2D_local = n_total_edges
      mesh%edge2D = edge2D_local
      mesh%edge2D_in = edge2D_in_local

      do n = 1, n_total_edges
        partit%myList_edge2D(n) = edge_perm(n)
      end do

      allocate(mesh%edges(2, edge2D_local))
      allocate(mesh%edge_tri(2, edge2D_local))

      do n = 1, n_total_edges
        eidx = edge_perm(n)
        mesh%edges(1, eidx) = edge_nodes(1, n)
        mesh%edges(2, eidx) = edge_nodes(2, n)
        mesh%edge_tri(1, eidx) = edge_elems(1, n)
        mesh%edge_tri(2, eidx) = edge_elems(2, n)
      end do

      deallocate(edge_nodes, edge_elems, edge_is_boundary, edge_perm)

      ! Re-sort myList_edge2D to identity after sorting
      do n = 1, edge2D_local
        partit%myList_edge2D(n) = n
      end do
    end if

    ! ========================================
    ! 1g. Set level arrays
    ! ========================================
    allocate(mesh%nlevels(elem2D_local))
    allocate(mesh%nlevels_nod2D(nod2D_total))
    allocate(mesh%ulevels(elem2D_local))
    allocate(mesh%ulevels_nod2D(nod2D_total))
    mesh%nlevels(:) = nl_in
    mesh%nlevels_nod2D(:) = nl_in
    mesh%ulevels(:) = 1
    mesh%ulevels_nod2D(:) = 1

    ! ========================================
    ! 1h. Build elem_edges
    ! ========================================
    allocate(mesh%elem_edges(3, elem2D_local))
    allocate(aux(elem2D_local))
    aux = 0
    do n = 1, edge2D_local
      do k = 1, 2
        q = mesh%edge_tri(k, n)
        if (q > 0 .and. q <= elem2D_local) then
          aux(q) = aux(q) + 1
          mesh%elem_edges(aux(q), q) = n
        end if
      end do
    end do
    deallocate(aux)

    ! Reorder elem_edges so that edge(q) does not contain node(q)
    do e = 1, elem2D_local
      elnodes = mesh%elem2D_nodes(:, e)
      eledges = mesh%elem_edges(:, e)
      do q = 1, 3
        do k = 1, 3
          if (mesh%edges(1, eledges(k)) /= elnodes(q) .and. &
              mesh%edges(2, eledges(k)) /= elnodes(q)) then
            mesh%elem_edges(q, e) = eledges(k)
            exit
          end if
        end do
      end do
    end do

    if (partit%mype == 0) then
      write(output_unit, '(A)') '  Mesh topology generated'
      write(output_unit, '(A,I8)') '  Internal edges: ', edge2D_in_local
      write(output_unit, '(A,I8)') '  Boundary edges: ', n_boundary
    end if

    ! ========================================
    ! 1i. Call existing mesh computation routines
    ! ========================================
    if (partit%mype == 0) write(output_unit, '(A)') '  Calling mesh computation routines...'

    call set_mesh_transform_matrix()
    call init_mpi_types(partit, mesh)
    call init_gatherLists(partit)
    call test_tri(partit, mesh)

    ! Build elem_neighbors and nod_in_elem2D directly instead of calling
    ! find_neighbors(), because the rectangular mesh has corner elements
    ! with only 1 neighbor (find_neighbors aborts if < 2).

    ! elem_neighbors: for each element, 3 neighbors from edge_tri
    allocate(mesh%elem_neighbors(3, elem2D_local))
    mesh%elem_neighbors = 0
    do elem = 1, elem2D_local
      eledges = mesh%elem_edges(:, elem)
      do j_nb = 1, 3
        elem1 = mesh%edge_tri(1, eledges(j_nb))
        if (elem1 == elem) elem1 = mesh%edge_tri(2, eledges(j_nb))
        mesh%elem_neighbors(j_nb, elem) = elem1
      end do
    end do

    ! nod_in_elem2D: for each node, list of elements containing it
    allocate(mesh%nod_in_elem2D_num(nod2D_total))
    mesh%nod_in_elem2D_num = 0
    do n = 1, elem2D_local
      do j_nb = 1, 3
        node = mesh%elem2D_nodes(j_nb, n)
        mesh%nod_in_elem2D_num(node) = mesh%nod_in_elem2D_num(node) + 1
      end do
    end do

    max_nod_in_elem = maxval(mesh%nod_in_elem2D_num)
    allocate(mesh%nod_in_elem2D(max_nod_in_elem, nod2D_total))
    mesh%nod_in_elem2D = 0
    mesh%nod_in_elem2D_num = 0
    do n = 1, elem2D_local
      do j_nb = 1, 3
        node = mesh%elem2D_nodes(j_nb, n)
        mesh%nod_in_elem2D_num(node) = mesh%nod_in_elem2D_num(node) + 1
        mesh%nod_in_elem2D(mesh%nod_in_elem2D_num(node), node) = n
      end do
    end do

    if (partit%mype == 0) write(output_unit, '(A)') '  Neighbor arrays built'

    call find_levels_min_e2n(partit, mesh)
    call mesh_areas(partit, mesh)
    call mesh_auxiliary_arrays(partit, mesh)

    ! ========================================
    ! 1j. Periodic post-remap
    ! ========================================
    if (is_periodic) then
      call periodic_post_remap()
    end if

    if (partit%mype == 0) then
      write(output_unit, '(A)') ''
      write(output_unit, '(A)') '================================================'
      write(output_unit, '(A)') 'Analytic mesh generation complete!'
      write(output_unit, '(A)') '================================================'
      write(output_unit, '(A)') ''
    end if

  contains

    !> Return node index for grid point (i,j) in the periodic mesh.
    !! i=1..nx, j=1..ny. Indices i=nx or j=ny produce halo node indices.
    function pnode(i, j) result(idx)
      integer, intent(in) :: i, j
      integer :: idx
      if (i <= nx-1 .and. j <= ny-1) then
        idx = (j-1)*(nx-1) + i
      else if (i == nx .and. j <= ny-1) then
        idx = nod2D_local + j              ! right halo
      else if (i <= nx-1 .and. j == ny) then
        idx = nod2D_local + (ny-1) + i     ! top halo
      else
        idx = nod2D_total                  ! corner halo
      end if
    end function pnode

    !> Post-remap for periodic mesh: fix geometry, merge areas, remap connectivity
    subroutine periodic_post_remap()
      integer :: h, r, nz
      integer, allocatable :: halo_to_real(:)
      real(kind=MP) :: half_Lx, half_Ly, val

      if (partit%mype == 0) write(output_unit, '(A)') '  Periodic post-remap...'

      ! Build halo-to-real mapping
      allocate(halo_to_real(eDim))
      ! Right halos: index nod2D_local+j -> real (j-1)*(nx-1)+1
      do h = 1, ny-1
        halo_to_real(h) = (h-1)*(nx-1) + 1
      end do
      ! Top halos: index nod2D_local+(ny-1)+i -> real i
      do h = 1, nx-1
        halo_to_real((ny-1) + h) = h
      end do
      ! Corner halo -> real 1
      halo_to_real(eDim) = 1

      ! Fix edge_cross_dxdy with minimum-image convention.
      ! For wrap-around edges, the element centroid is on the opposite side of
      ! the domain from the edge center. The raw cross-dxdy points across the
      ! whole domain instead of wrapping. Correct using minimum-image.
      half_Lx = 0.5_MP * Lx
      half_Ly = 0.5_MP * Ly
      do n = 1, mesh%edge2D
        do k = 1, 3, 2  ! slots 1,3 = x-component; k+1 = y-component
          val = mesh%edge_cross_dxdy(k, n)
          if (abs(val) > half_Lx) then
            mesh%edge_cross_dxdy(k, n) = val - sign(Lx, val)
          end if
          val = mesh%edge_cross_dxdy(k+1, n)
          if (abs(val) > half_Ly) then
            mesh%edge_cross_dxdy(k+1, n) = val - sign(Ly, val)
          end if
        end do
      end do

      ! Merge halo node areas into real nodes
      do h = 1, eDim
        r = halo_to_real(h)
        do nz = 1, mesh%nl - 1
          mesh%area(nz, r) = mesh%area(nz, r) + mesh%area(nz, nod2D_local + h)
          mesh%areasvol(nz, r) = mesh%areasvol(nz, r) + mesh%areasvol(nz, nod2D_local + h)
        end do
      end do

      ! Recompute inverses and mesh resolution for real nodes
      do n = 1, nod2D_local
        do nz = 1, mesh%nl - 1
          if (mesh%area(nz, n) > 0.0_MP) then
            mesh%area_inv(nz, n) = 1.0_MP / mesh%area(nz, n)
          else
            mesh%area_inv(nz, n) = 0.0_MP
          end if
          if (mesh%areasvol(nz, n) > 0.0_MP) then
            mesh%areasvol_inv(nz, n) = 1.0_MP / mesh%areasvol(nz, n)
          else
            mesh%areasvol_inv(nz, n) = 0.0_MP
          end if
        end do
        mesh%mesh_resolution(n) = sqrt(mesh%areasvol(mesh%ulevels_nod2D(n), n) / pi) * 2.0_MP
      end do

      ! Remap elem2D_nodes: replace halo indices with real indices
      do n = 1, mesh%elem2D
        do k = 1, 3
          h = mesh%elem2D_nodes(k, n)
          if (h > nod2D_local) then
            mesh%elem2D_nodes(k, n) = halo_to_real(h - nod2D_local)
          end if
        end do
      end do

      ! Remap edges: replace halo indices with real indices
      do n = 1, mesh%edge2D
        do k = 1, 2
          h = mesh%edges(k, n)
          if (h > nod2D_local) then
            mesh%edges(k, n) = halo_to_real(h - nod2D_local)
          end if
        end do
      end do

      ! Set eDim_nod2D = 0 — subsequent allocations use only real nodes
      partit%eDim_nod2D = 0

      deallocate(halo_to_real)

      if (partit%mype == 0) write(output_unit, '(A)') '  Periodic post-remap complete'
    end subroutine periodic_post_remap

  end subroutine generate_analytic_mesh

end module analytic_mesh_module
