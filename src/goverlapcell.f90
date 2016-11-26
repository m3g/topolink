!
! Subroutine that computes the gradient of the overlap in a cell
!

subroutine goverlapcell(n,i,ix,iy,iz,x,ibox,jbox,kbox,g)

  use functionpars
  use linkedcells
  implicit none
  integer :: n, i, ix, iy, iz, ibox, jbox, kbox, iatom
  double precision :: x(n), g(n), dx, dy, dz, d, dconst, vdwrad2

  if ( ibox < 1 .or. jbox < 1 .or. kbox < 1 ) return
  if ( ibox > nboxesx .or. jbox > nboxesy .or. kbox > nboxesz ) return
  
  vdwrad2 = sigma(i)**2
  iatom = ifirstbox(ibox,jbox,kbox) 
  do while( iatom /= 0 )
    
    if ( skip(iatom) ) then
      iatom = inextbox(iatom)
      cycle
    end if

    dx = x(ix)-coor(iatom,1)
    dy = x(iy)-coor(iatom,2)
    dz = x(iz)-coor(iatom,3)
    d = dx**2 + dy**2 + dz**2

    if ( d < vdwrad2 ) then
      d = dsqrt(d)
      dconst = 2.d0*kvdw*(sigma(i)-d)/d
      g(ix) = g(ix) - dconst*dx
      g(iy) = g(iy) - dconst*dy
      g(iz) = g(iz) - dconst*dz
    end if

    iatom = inextbox(iatom)
  end do

end subroutine goverlapcell
