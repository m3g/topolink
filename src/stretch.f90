!
! Subroutine that computes the stretch energy of the linker
!

double precision function stretch(n,x)

  use functionpars
  implicit none
  integer :: n, i, ix, iy, iz, jx, jy, jz
  double precision :: x(n), d, d2

  stretch = 0.d0

  ! First bond, to fixed atom1

  d2 = (x(1)-coor(atom1,1))**2 + &
       (x(2)-coor(atom1,2))**2 + &
       (x(3)-coor(atom1,3))**2
  if ( d2 <= dbond2 ) then
    stretch = stretch + d2
  else
    d = dsqrt(d2)
    stretch = stretch + d2 + kbond*(d-dbond)**2
  end if

  ! Bonds within the linker, with variable positions for both atoms

  do i = 1, nlinkatoms-1
    ix = (i-1)*3 + 1
    iy = ix + 1
    iz = ix + 2
    jx = ix + 3
    jy = jx + 1
    jz = jx + 2
    d2 = (x(jx)-x(ix))**2 + &
         (x(jy)-x(iy))**2 + &
         (x(jz)-x(iz))**2
    if ( d2 <= dbond2 ) then
      stretch = stretch + d2
    else
      d = dsqrt(d2)
      stretch = stretch + d2 + kbond*(d-dbond)**2
    end if
  end do

  ! Last bond, to fixed atom2 

  ix = (nlinkatoms-1)*3 + 1
  iy = ix + 1
  iz = ix + 2
  d2 = (coor(atom2,1)-x(ix))**2 + &
       (coor(atom2,2)-x(iy))**2 + &
       (coor(atom2,3)-x(iz))**2
  if ( d2 <= dbond2 ) then
    stretch = stretch + d2
  else
    d = dsqrt(d2)
    stretch = stretch + d2 + kbond*(d-dbond)**2
  end if

end function stretch
