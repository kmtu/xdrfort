!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!

program example

    ! 1. Use the xdr interface
    use xdr, only: xtcfile

    implicit none

    ! 2. Declare a variable of type xtcfile
    type(xtcfile) :: xtcf
    type(xtcfile) :: xtc_out

    ! 3. Initialize it with the names of xtc files you want to read in and write out
    call xtcf % init("example.xtc")
    call xtc_out % init("example_out.xtc", 'w')

    ! 4. Read in each configuration. Everything is stored in the xtcfile type (precision, time,
    !    step, no of atoms, positions, etc.). Look in the xtc module for more details.
    !    You can save the positions in the loop for your calculations in another array, or 
    !    do your calculations after each read.

    call xtcf % read

    do while ( xtcf % STAT == 0 )

        ! Just an example to show what was read in
        write(*,'(a,f12.6,a,i0)') " Time (ps): ", xtcf % time, "  Step: ", xtcf % STEP
        write(*,'(a,f12.6,a,i0)') " Precision: ", xtcf % prec, "  No. Atoms: ", xtcf % NATOMS
        write(*,'(3f9.3)') xtcf % pos

        ! This is the same order as found in the GRO format fyi
        write(*,'(11f9.5)') xtcf % box(1,1), xtcf % box(2,2), xtcf % box(3,3), &
                            xtcf % box(1,2), xtcf % box(1,3), & 
                            xtcf % box(2,1), xtcf % box(2,3), &
                            xtcf % box(3,1), xtcf % box(3,2) 
        call xtc_out % write(xtcf % natoms, xtcf % step, xtcf % time, xtcf % box, xtcf % pos, xtcf % prec)

        call xtcf % read
    end do

    ! 5. Close the file
    call xtcf % close
    call xtc_out % close

end program example
