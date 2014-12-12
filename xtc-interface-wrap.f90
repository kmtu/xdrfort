!  XDR Fortran Interface with Wrappers
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/

module xtc

    use, intrinsic :: iso_c_binding, only: C_PTR, C_CHAR, C_FLOAT, C_INT

    implicit none
    private

!   *** xtcfile type
!   box     - triclinic pbc box of the configuration
!   NATOMS  - number of atoms in the configuration.
!   pos     - positions read in (3,NATOMS)
!   prec    - precision of the coordinates read in
!   STEP    - step number of configuration.
!   STAT    - status of operation. 0 = good
!   time    - time of the configuration
!   xd      - pointer from libxdrfile.

!   Should always call init first. Then call read in a loop and do your
!   calculations. After the loops call close.

    type, public :: xtcfile
      type(xdrfile), pointer :: xd
      real(C_FLOAT), allocatable :: pos(:,:)
      integer(C_INT) :: NATOMS, STEP, STAT
      real(C_FLOAT) :: box(3,3), prec, time
    contains
      procedure :: init => init_xtc
      procedure :: read => read_xtc_wrap
      procedure :: close => close_xtc
    end type

    ! the data type located in libxdrfile
    type, bind(C) :: xdrfile
      type(C_PTR) :: fp, xdr
      character(kind=C_CHAR) :: mode
      integer(C_INT) :: buf1, buf1size, buf2, buf2size
    end type xdrfile

    ! interface with libxdrfile
    interface 

      integer(C_INT) function read_xtc_natoms(filename,NATOMS) bind(C, name='read_xtc_natoms')
        import
        character(kind=C_CHAR), intent(in) :: filename
        integer(C_INT), intent(out) :: NATOMS
      end function

      type(C_PTR) function xdrfile_open(filename,mode) bind(C, name='xdrfile_open')
        import
        character(kind=C_CHAR), intent(in) :: filename(*), mode(*)
      end function

      integer(C_INT) function read_xtc(xd,NATOMS,STEP,time,box,x,prec) bind(C, name='read_xtc')
        import
        type(xdrfile), intent(in) :: xd
        integer(C_INT), intent(out) :: NATOMS, STEP
        real(C_FLOAT), intent(out) :: time, prec, box(*), x(*)
      end function

      integer(C_INT) function xdrfile_close(xd) bind(C,name='xdrfile_close')
        import
        type(xdrfile), intent(in) :: xd
      end function

    end interface

contains

    ! our wrappers for the xtcfile class
    subroutine init_xtc(xtc,filename_in)

        use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_CHAR, c_f_pointer

        implicit none
        class(xtcfile), intent(inout) :: xtc
        type(C_PTR) :: xd_c
        character (len=*), intent(in) :: filename_in
        character (len=206) :: filename
        logical :: ex

        inquire(file=trim(filename_in),exist=ex)

        if (ex .eqv. .false.) then
            write(0,*)
            write(0,'(a)') " Error: "//trim(filename_in)//" does not exist."
            write(0,*)
            stop
        end if

        ! Set the file name to be read in for C.
        filename = trim(filename_in)//C_NULL_CHAR

        ! Get number of atoms in system and allocate position array.
        xtc % STAT = read_xtc_natoms(filename,xtc % NATOMS)
        allocate(xtc % pos(3,xtc % NATOMS))

        if (xtc % STAT /= 0) then
            write(0,*)
            write(0,'(a)') " Error reading in "//trim(filename_in)//". Is it really an xtc file?"
            write(0,*)
            stop
        end if

        ! Open the file for reading. Convert C pointer to Fortran pointer.
        xd_c = xdrfile_open(filename,"r")
        call c_f_pointer(xd_c,xtc % xd)

        write(0,'(a)') " Opened "//trim(filename)//" for reading."
        write(0,'(a,i0,a)') " ",xtc % NATOMS, " atoms present in system."
        write(0,*)

    end subroutine init_xtc

    subroutine read_xtc_wrap(xtc)

        implicit none
        class(xtcfile), intent(inout) :: xtc
        real :: box_trans(3,3)

        xtc % STAT = read_xtc(xtc % xd,xtc % NATOMS,xtc % STEP,xtc % time,box_trans,xtc % pos,xtc % prec)

        ! C is row-major, whereas Fortran is column major. Hence the following.
        xtc % box = transpose(box_trans)

    end subroutine read_xtc_wrap

    subroutine close_xtc(xtc)

        implicit none
        class(xtcfile), intent(inout) :: xtc

        xtc % STAT = xdrfile_close(xtc % xd)
        deallocate(xtc % pos)

    end subroutine close_xtc

end module xtc
