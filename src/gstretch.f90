!
! Subroutine that adds to the gradient the component associated with
! the stretch of the linker
!

subroutine gstretch(n,x,g)

  use functionpars
  implicit none
  integer :: n, i, ix, iy, iz, jx, jy, jz
  double precision :: x(*), g(*), d2, dx, dy, dz, d, gx, gy, gz

  ! First bond, to fixed atom1

  dx = x(1) - coor(atom1,1)
  dy = x(2) - coor(atom1,2)
  dz = x(3) - coor(atom1,3)
  d2 = dx**2 + dy**2 + dz**2
  if ( d2 <= dbond2 ) then
    g(1) = g(1) + 2.d0*dx
    g(2) = g(2) + 2.d0*dy
    g(3) = g(3) + 2.d0*dz
  else
    d = dsqrt(d2)
    g(1) = g(1) + 2.d0*( dx + kbond*(d-dbond)*dx/d )
    g(2) = g(2) + 2.d0*( dy + kbond*(d-dbond)*dy/d )
    g(3) = g(3) + 2.d0*( dz + kbond*(d-dbond)*dz/d )
  end if

  ! Normal bonds (not first neither last)
  
  do i = 1, nlinkatoms-1

    ix = (i-1)*3 + 1
    iy = ix + 1
    iz = ix + 2
    jx = ix + 3
    jy = jx + 1
    jz = jx + 2
    dx = x(jx)-x(ix)
    dy = x(jy)-x(iy)
    dz = x(jz)-x(iz)
    d2 = dx**2 + dy**2 + dz**2

    if ( d2 <= dbond2 ) then
     gx = 2.d0*dx
     gy = 2.d0*dy
     gz = 2.d0*dz
    else
      d = dsqrt(d2)
      gx = 2.d0*( dx + kbond*(d-dbond)*dx/d )
      gy = 2.d0*( dy + kbond*(d-dbond)*dy/d )
      gz = 2.d0*( dz + kbond*(d-dbond)*dz/d )
    end if
   
    g(ix) = g(ix) - gx
    g(iy) = g(iy) - gy
    g(iz) = g(iz) - gz

    g(jx) = g(jx) + gx
    g(jy) = g(jy) + gy
    g(jz) = g(jz) + gz

  end do

  ! Last bond, to fixed atom2

  ix = (nlinkatoms-1)*3 + 1
  iy = ix + 1
  iz = ix + 2
  dx = coor(atom2,1) - x(ix)
  dy = coor(atom2,2) - x(iy)
  dz = coor(atom2,3) - x(iz)
  d2 = dx**2 + dy**2 + dz**2
  if ( d2 <= dbond2 ) then
    g(ix) = g(ix) - 2.d0*dx
    g(iy) = g(iy) - 2.d0*dy
    g(iz) = g(iz) - 2.d0*dz
  else
    d = dsqrt(d2)
    g(ix) = g(ix) - 2.d0*( dx + kbond*(d-dbond)*dx/d )
    g(iy) = g(iy) - 2.d0*( dy + kbond*(d-dbond)*dy/d )
    g(iz) = g(iz) - 2.d0*( dz + kbond*(d-dbond)*dz/d )
  end if

end subroutine gstretch
