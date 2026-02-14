!> @file mesh_output.F90
!! @brief Write mesh arrays as raw binary files for Python visualization
!! @details Writes coord_nod2D, elem2D_nodes, zbar, and depth as raw binary
!!          files (access='stream') that numpy.fromfile() can read directly.
!!          A small mesh_info.txt file provides array dimensions.

module mesh_output_module
  use MOD_MESH
  use MOD_TRACER
  use o_PARAM, only: WP, MP, r_earth
  use fortran_utils, only: mkdir
  implicit none
  private
  public :: write_mesh_for_python
  public :: open_scalar_output, write_scalar_step, close_scalar_output

contains

  subroutine write_mesh_for_python(mesh, output_dir)
    type(t_mesh), intent(in)   :: mesh
    character(*), intent(in)   :: output_dir

    integer :: iunit, i, n
    real(8), allocatable    :: coords(:,:), zbar_f64(:), depth_f64(:)
    integer(4), allocatable :: elems(:,:)

    ! Create output directory (ignores error if already exists)
    call mkdir(output_dir)

    ! ---- mesh_info.txt ----
    open(newunit=iunit, file=output_dir//'/mesh_info.txt', &
         status='replace', action='write')
    write(iunit, '(A,I0)') 'nod2D=', mesh%nod2D
    write(iunit, '(A,I0)') 'elem2D=', mesh%elem2D
    write(iunit, '(A,I0)') 'nl=', mesh%nl
    close(iunit)

    ! ---- coord_nod2D.bin : float64 (2, nod2D) column-major, meters ----
    allocate(coords(2, mesh%nod2D))
    do n = 1, mesh%nod2D
      coords(1, n) = real(mesh%coord_nod2D(1, n), 8) * real(r_earth, 8)
      coords(2, n) = real(mesh%coord_nod2D(2, n), 8) * real(r_earth, 8)
    end do
    open(newunit=iunit, file=output_dir//'/coord_nod2D.bin', &
         access='stream', status='replace', action='write', form='unformatted')
    write(iunit) coords
    close(iunit)
    deallocate(coords)

    ! ---- elem2D_nodes.bin : int32 (3, elem2D) column-major, 0-based ----
    allocate(elems(3, mesh%elem2D))
    do n = 1, mesh%elem2D
      do i = 1, 3
        elems(i, n) = int(mesh%elem2D_nodes(i, n) - 1, 4)
      end do
    end do
    open(newunit=iunit, file=output_dir//'/elem2D_nodes.bin', &
         access='stream', status='replace', action='write', form='unformatted')
    write(iunit) elems
    close(iunit)
    deallocate(elems)

    ! ---- zbar.bin : float64 (nl,) ----
    allocate(zbar_f64(mesh%nl))
    do i = 1, mesh%nl
      zbar_f64(i) = real(mesh%zbar(i), 8)
    end do
    open(newunit=iunit, file=output_dir//'/zbar.bin', &
         access='stream', status='replace', action='write', form='unformatted')
    write(iunit) zbar_f64
    close(iunit)
    deallocate(zbar_f64)

    ! ---- depth.bin : float64 (nod2D,) ----
    allocate(depth_f64(mesh%nod2D))
    do n = 1, mesh%nod2D
      depth_f64(n) = real(mesh%depth(n), 8)
    end do
    open(newunit=iunit, file=output_dir//'/depth.bin', &
         access='stream', status='replace', action='write', form='unformatted')
    write(iunit) depth_f64
    close(iunit)
    deallocate(depth_f64)

    write(*, '(A,A)') '  Mesh written to ', output_dir

  end subroutine write_mesh_for_python

  subroutine open_scalar_output(mesh, nod2D, nz, nsteps, num_tracers, output_dir, file_units)
    type(t_mesh), intent(in)   :: mesh
    integer, intent(in)        :: nod2D, nz, nsteps, num_tracers
    character(*), intent(in)   :: output_dir
    integer, intent(out)       :: file_units(num_tracers)

    integer :: iunit, n
    character(len=32) :: tracer_names(2)

    tracer_names(1) = 'temperature'
    tracer_names(2) = 'salinity'

    ! Create output directory (ignores error if already exists)
    call mkdir(output_dir)

    ! Write scalar_info.txt
    open(newunit=iunit, file=output_dir//'/scalar_info.txt', &
         status='replace', action='write')
    write(iunit, '(A,I0)') 'nod2D=', nod2D
    write(iunit, '(A,I0)') 'nz=', nz
    write(iunit, '(A,I0)') 'nsteps=', nsteps
    write(iunit, '(A,I0)') 'num_tracers=', num_tracers
    write(iunit, '(A,I0)') 'wp=', WP
    do n = 1, num_tracers
      write(iunit, '(A,I0,A,A)') 'tracer_', n, '=', trim(tracer_names(n))
    end do
    close(iunit)

    ! Open one stream file per tracer
    do n = 1, num_tracers
      open(newunit=file_units(n), &
           file=output_dir//'/'//trim(tracer_names(n))//'.bin', &
           access='stream', status='replace', action='write', form='unformatted')
    end do

    write(*, '(A,A)') '  Scalar output files opened in ', output_dir

  end subroutine open_scalar_output


  subroutine write_scalar_step(file_units, num_tracers, tracers, nz, nod2D)
    integer, intent(in)        :: file_units(:)
    integer, intent(in)        :: num_tracers, nz, nod2D
    type(t_tracer), intent(in) :: tracers

    integer :: n
    real(8), allocatable :: buf(:,:)

    allocate(buf(nz, nod2D))

    do n = 1, num_tracers
      buf(:,:) = real(tracers%data(n)%values(1:nz, 1:nod2D), 8)
      write(file_units(n)) buf
    end do

    deallocate(buf)

  end subroutine write_scalar_step


  subroutine close_scalar_output(file_units, num_tracers)
    integer, intent(in) :: file_units(:)
    integer, intent(in) :: num_tracers

    integer :: n

    do n = 1, num_tracers
      close(file_units(n))
    end do

    write(*, '(A)') '  Scalar output files closed'

  end subroutine close_scalar_output

end module mesh_output_module
