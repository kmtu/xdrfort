!  XDR Fortran Interface with Wrappers
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/

!  Modified by Kai-Min Tu (2014)
!  Based on Barnett's work, the interface for TRR file is incorporated
!  https://github.com/kmtu/

module xdr

    use, intrinsic :: iso_c_binding, only: C_PTR, C_CHAR, C_FLOAT, C_INT

    implicit none
    private

    type, abstract :: trjfile
      type(xdrfile), pointer :: xd
      integer(C_INT) :: natoms, step, stat
      real(C_FLOAT) :: box(3,3), time
      character(len=1) :: mode
    contains
      procedure :: init => init_xdr
      procedure :: read => read_xdr
      procedure :: close => close_xdr
    end type

!   *** xtcfile type
!   box     - triclinic pbc box of the configuration
!   natoms  - number of atoms in the configuration.
!   pos     - positions read in (3,natoms)
!   prec    - precision of the coordinates read in
!   step    - step number of configuration.
!   stat    - status of operation. 0 = good
!   time    - time of the configuration
!   xd      - pointer from libxdrfile.

!   Should always call init first. Then call read in a loop and do your
!   calculations. After the loops call close.

    type, extends(trjfile), public :: xtcfile
      real(C_FLOAT), allocatable :: pos(:,:)
      real(C_FLOAT) :: prec
    contains
      procedure :: write => write_xtcfile
    end type

!   *** trrfile type
!   box     - triclinic pbc box of the configuration
!   natoms  - number of atoms in the configuration.
!   pos     - positions read in (3,natoms)
!   vel     - velocities read in (3,natoms)
!   force   - forces read in (3,natoms)
!   lambda  - lambda value for free energy perturbation calculations
!   step    - step number of configuration.
!   stat    - status of operation. 0 = good
!   time    - time of the configuration
!   xd      - pointer from libxdrfile.

    type, extends(trjfile), public :: trrfile
      real(C_FLOAT), allocatable :: pos(:,:), vel(:,:), force(:,:)
      real(C_FLOAT) :: lambda
    contains
      procedure :: write => write_trrfile
    end type

    ! the data type located in libxdrfile
    type, bind(C) :: xdrfile
      type(C_PTR) :: fp, xdr
      character(kind=C_CHAR) :: mode
      type(C_PTR) :: buf1, buf2
      integer(C_INT) :: buf1size, buf2size
    end type xdrfile

    ! interface with libxdrfile
    interface 

      type(C_PTR) function xdrfile_open(filename,mode) bind(C, name='xdrfile_open')
        import
        character(kind=C_CHAR), intent(in) :: filename(*), mode(*)
      end function

      integer(C_INT) function xdrfile_close(xd) bind(C,name='xdrfile_close')
        import
        type(xdrfile), intent(in) :: xd
      end function

      ! xtc
      integer(C_INT) function read_xtc_natoms(filename,natoms) bind(C)
        import
        character(kind=C_CHAR), intent(in) :: filename
        integer(C_INT), intent(out) :: natoms
      end function

      integer(C_INT) function read_xtc(xd,natoms,step,time,box,x,prec) bind(C)
        import
        type(xdrfile), intent(in) :: xd
        integer(C_INT), intent(in), value :: natoms
        integer(C_INT), intent(out) :: step
        real(C_FLOAT), intent(out) :: time, prec, box(*), x(*)
      end function

      integer(C_INT) function write_xtc(xd,natoms,step,time,box,x,prec) bind(C)
        import
        type(xdrfile), intent(in) :: xd
        integer(C_INT), intent(in), value :: natoms, step
        real(C_FLOAT), intent(in), value :: time, prec
        real(C_FLOAT), intent(in) :: box(*), x(*)
      end function

      ! trr
      integer(C_INT) function read_trr_natoms(filename,natoms) bind(C)
        import
        character(kind=C_CHAR), intent(in) :: filename
        integer(C_INT), intent(out) :: natoms
      end function

      integer(C_INT) function read_trr(xd,natoms,step,time,lambda,box,x,v,f) bind(C)
        import
        type(xdrfile), intent(in) :: xd
        integer(C_INT), intent(in), value :: natoms
        integer(C_INT), intent(out) :: step
        real(C_FLOAT), intent(out) :: time, lambda, box(*), x(*), v(*), f(*)
      end function

      integer(C_INT) function write_trr(xd,natoms,step,time,lambda,box,x,v,f) bind(C)
        import
        type(xdrfile), intent(in) :: xd
        integer(C_INT), intent(in), value :: natoms, step
        real(C_FLOAT), intent(in), value :: time, lambda
        real(C_FLOAT), intent(in) :: box(*), x(*), v(*), f(*)
      end function

    end interface

contains

    ! our wrappers for the trjfile class
    subroutine init_xdr(trj,filename_in,mode_opt)

        use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_CHAR, c_f_pointer

        implicit none
        class(trjfile), intent(inout) :: trj
        type(C_PTR) :: xd_c
        character (len=*), intent(in) :: filename_in
        character (len=206) :: filename
        logical :: ex
        character(len=1), optional, intent(in) :: mode_opt

        if (present(mode_opt)) then
          trj % mode = mode_opt
        else
          trj % mode = 'r'
        end if

        ! Set the file name to be read in for C.
        filename = trim(filename_in)//C_NULL_CHAR

        select case (trj % mode)
        case ('r')
          inquire(file=trim(filename_in),exist=ex)

          if (ex .eqv. .false.) then
              write(0,*)
              write(0,'(a)') " Error: "//trim(filename_in)//" does not exist."
              write(0,*)
              stop
          end if

          select type (trj)
          type is (xtcfile)
            ! Get number of atoms in system and allocate position array.
            trj % stat = read_xtc_natoms(filename,trj % natoms)

            if (trj % stat /= 0) then
                write(0,*)
                write(0,'(a)') " Error reading in "//trim(filename_in)//". Is it really an xtc file?"
                write(0,*)
                stop
            end if

            allocate(trj % pos(3,trj % natoms))

          type is (trrfile)
            ! Get number of atoms in system and allocate position, velocity, force arrays.
            trj % stat = read_trr_natoms(filename,trj % natoms)

            if (trj % stat /= 0) then
                write(0,*)
                write(0,'(a)') " Error reading in "//trim(filename_in)//". Is it really an trr file?"
                write(0,*)
                stop
            end if

            allocate(trj % pos(3,trj % natoms))
            allocate(trj % vel(3,trj % natoms))
            allocate(trj % force(3,trj % natoms))
          end select

          write(0,'(a)') " Open "//trim(filename)//" for reading."
          write(0,'(a,i0,a)') " ",trj % natoms, " atoms present in system."
          write(0,*)

        case ('w')
          write(0,'(a)') " Open "//trim(filename)//" for writing."
          write(0,*)

        case default
          write(0,*)
          write(0,*) "Unknown file mode: '"//trj % mode//"'. It can only be 'r' or 'w'."
          write(0,*)
          stop
        end select

        ! Open the file for reading or writing. Convert C pointer to Fortran pointer.
        xd_c = xdrfile_open(filename, trj % mode)
        call c_f_pointer(xd_c,trj % xd)
    end subroutine init_xdr

    subroutine read_xdr(trj)

        implicit none
        class(trjfile), intent(inout) :: trj

        select type (trj)

        type is (xtcfile)
          trj % stat = read_xtc(trj % xd,trj % natoms,trj % step,trj % time,trj % box,trj % pos,trj % prec)

        type is (trrfile)
          trj % stat = read_trr(trj % xd,trj % natoms,trj % step,trj % time,trj % lambda, trj % box,trj % pos,trj % vel,trj % force)

        end select

    end subroutine read_xdr

    subroutine write_xtcfile(xtc, natoms, step, time, box, pos, prec)
        implicit none
        class(xtcfile), intent(inout) :: xtc
        integer, intent(in) :: natoms, step
        real, intent(in) :: time, box(3,3), pos(:, :), prec

        xtc % stat = write_xtc(xtc % xd, natoms, step, time, box, pos, prec)
    end subroutine write_xtcfile

    subroutine write_trrfile(trr, natoms, step, time, lambda, box, pos, vel, force)
        implicit none
        class(trrfile), intent(inout) :: trr
        integer, intent(in) :: natoms, step
        real, intent(in) :: time, lambda, box(3,3), pos(*), vel(*), force(*)

        trr % stat = write_trr(trr % xd, natoms, step, time, lambda, box, pos, vel, force)
    end subroutine write_trrfile

    subroutine close_xdr(trj)

        implicit none
        class(trjfile), intent(inout) :: trj

        if (trj % mode == 'r') then
          select type (trj)
          type is (xtcfile)
            deallocate(trj % pos)

          type is (trrfile)
            deallocate(trj % pos)
            deallocate(trj % vel)
            deallocate(trj % force)
          end select
        end if 

        trj % stat = xdrfile_close(trj % xd)
    end subroutine close_xdr

end module xdr
