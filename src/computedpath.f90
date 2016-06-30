!
! Function that computes the path distance
!

double precision function computedpath(n,x)

  use functionpars
  implicit none
  integer :: n, i, ix, iy, iz, jx, jy, jz
  double precision :: x(n)


  computedpath = dsqrt( (x(1)-coor(atom1,1))**2 + &
                        (x(2)-coor(atom1,2))**2 + &
                        (x(3)-coor(atom1,3))**2 )
  do i = 1, nlinkatoms-1
    ix = (i-1)*3+1
    iy = ix + 1
    iz = ix + 2
    jx = ix + 3
    jy = jx + 1
    jz = jx + 2
    computedpath = computedpath + dsqrt( (x(jx)-x(ix))**2 + &
                                         (x(jy)-x(iy))**2 + &
                                         (x(jz)-x(iz))**2 )
  end do
  ix = (nlinkatoms-1)*3 + 1
  iy = ix + 1
  iz = ix + 2
  computedpath = computedpath + dsqrt( (coor(atom2,1)-x(ix))**2 + &
                                       (coor(atom2,2)-x(iy))**2 + &
                                       (coor(atom2,3)-x(iz))**2 )

end function computedpath
