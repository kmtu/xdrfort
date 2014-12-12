!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!

program read_xtc_prog

    ! 1. Use the xdr interface
    use xtc, only: xtcfile

    implicit none

    ! 2. Declare a variable of type xtcfile
    type(xtcfile) :: xtc

    ! 3. Initialize it with the name of xtc file you want to read in.
    call xtc % init("traj.xtc")

    ! 4. Read in each configuration. Everything is stored in the xtcfile type (precision, time,
    !    step, no of atoms, positions, etc.). Look in the xtc module for more details.
    !    You can save the positions in the loop for your calculations in another array, or 
    !    do your calculations after each read.

    call xtc % read

    do while ( xtc % STAT == 0 )

        ! Just an example to show what was read in
        write(*,'(a,f12.6,a,i0)') " Time (ps): ", xtc % time, "  Step: ", xtc % STEP
        write(*,'(a,f12.6,a,i0)') " Precision: ", xtc % prec, "  No. Atoms: ", xtc % NATOMS
        write(*,'(3f9.3)') xtc % pos

        ! This is the same order as found in the GRO format fyi
        write(*,'(11f9.5)') xtc % box(1,1), xtc % box(2,2), xtc % box(3,3), &
                            xtc % box(1,2), xtc % box(1,3), & 
                            xtc % box(2,1), xtc % box(2,3), &
                            xtc % box(3,1), xtc % box(3,2) 
        call xtc % read

    end do

    ! 5. Close the file
    call xtc % close

end program read_xtc_prog
