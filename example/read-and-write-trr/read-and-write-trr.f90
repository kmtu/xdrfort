!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!
!  Modified by Kai-Min Tu (2014)
!  Based on Barnett's work, the interface for TRR file is incorporated
!  https://github.com/kmtu/

program example

    ! 1. Use the xdr interface
    use xdr, only: trrfile

    implicit none

    ! 2. Declare a variable of type trrfile
    type(trrfile) :: trr
    type(trrfile) :: trr_out

    ! 3. Initialize it with the names of trr files you want to read in and write out
    call trr % init("example.trr")
    call trr_out % init("example_out.trr", 'w')

    ! 4. Read in each configuration. Everything is stored in the trrfile type (lambda, time,
    !    step, no of atoms, positions, etc.). Look in the xtc module for more details.
    !    You can save the positions in the loop for your calculations in another array, or 
    !    do your calculations after each read.

    call trr % read

    do while ( trr % STAT == 0 )
        ! Just an example to show what was read in
        write(*,'(a,f12.6,a,i0)') " Time (ps): ", trr % time, "  Step: ", trr % STEP
        write(*,'(a,f12.6,a,i0)') " Lambda: ", trr % lambda, "  No. Atoms: ", trr % NATOMS
        write(*,*)

        write(*,*) "=== box ==="
        ! This is the same order as found in the GRO format fyi
        write(*,'(11f9.5)') trr % box(1,1), trr % box(2,2), trr % box(3,3), &
                            trr % box(1,2), trr % box(1,3), & 
                            trr % box(2,1), trr % box(2,3), &
                            trr % box(3,1), trr % box(3,2) 
        write(*,*)

        write(*,*) "=== positions ==="
        write(*,'(3f9.3)') trr % pos
        write(*,*)
        write(*,*) "=== velocities ==="
        write(*,'(3f9.3)') trr % vel
        write(*,*)
        write(*,*) "=== forces ==="
        write(*,'(3f9.3)') trr % force
        write(*,*)

        call trr_out % write(trr % natoms, trr % step, trr % time, trr % lambda, trr % box, trr % pos, trr % vel, trr % force)
        if (trr_out % stat /= 0) then
          write(0,*) "Error: something was wrong when writing data"
        end if

        call trr % read
    end do

    ! 5. Close the file
    call trr % close
    call trr_out % close

end program example
