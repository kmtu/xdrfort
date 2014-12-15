! This program calculates the ow-ow rdf
! for pure water in NVT ensemble with periodic cubic boundaries
program rdf
    use xdr, only: xtcfile
    implicit none
    type(xtcfile) :: xtc
    real(kind=8), parameter :: pi = 3.141592653589793238462643383
    real(kind=8), parameter :: dr = 0.002 !nm
    real(kind=8), allocatable :: pos(:,:), g(:), r(:)
    real(kind=8) :: xr(3), dxr  !dxr: distance between two atoms
    real(kind=8) :: r_max, box_dim, dv, rho
    integer :: nhist, ng, ig, npos, i, j

    call xtc % init("w100-nvt.xtc")

    npos = xtc % NATOMS / 3  !number of water molecules (NATOMS is obtained after calling init)
    allocate(pos(3, npos))

    call xtc % read
    
    ! box information cannot be obtained until at least one read call
    box_dim = xtc % box(1,1)
    r_max = box_dim / 2d0
    nhist = ceiling(r_max / dr)
    allocate(g(nhist))
    allocate(r(nhist))
    g = 0d0
    ng = 0
    r = [((i - 0.5) * dr, i = 1, nhist)]  !set r-scales as the middle points (Fortran comprehension list)

    do while ( xtc % STAT == 0 )
        ng = ng + 1
        pos = xtc % pos(:, 1:xtc % NATOMS:3)  !get the position of OW (every 3rd atom)
        do i = 1, npos
          do j = 1, npos
            if (i /= j) then
              xr = pos(:, i) - pos(:, j)
              xr = abs(xr - box_dim * nint(xr / box_dim))  !wrap distance (periodic boundary condition)
              dxr = sqrt(sum(xr**2))
              ig = ceiling(dxr / dr)
              if (ig <= nhist) then
                if (ig == 0) then
                  ig = 1  !index of g(:) begins from 1
                end if
                g(ig) = g(ig) + 1
              end if
            end if
          end do
        end do
        call xtc % read
    end do

    ! 5. Close the file
    call xtc % close

    ! normalize rdf
    rho = dble(npos) / box_dim**3  !number density
    do i = 1, nhist
      dv = (4d0 / 3d0) * pi * (i**3 - (i-1)**3) * dr**3
      g(i) = g(i) / (ng * npos * dv * rho)
    end do

    ! output results
    do i = 1, nhist
      write(*,'(2f7.3)') r(i), g(i)
    end do
end program rdf
