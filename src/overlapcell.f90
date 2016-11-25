!
! Subroutine that computes the overlap in a cell
!

subroutine overlapcell(x,y,z,ibox,jbox,kbox,overlap)

  use functionpars
  use linkedcells
  implicit none
  integer :: ibox, jbox, kbox, iatom
  double precision :: x, y, z, overlap, d

  if ( ibox < 1 .or. jbox < 1 .or. kbox < 1 ) return
  if ( ibox > nboxesx .or. jbox > nboxesy .or. kbox > nboxesz ) return
  
  iatom = ifirstbox(ibox,jbox,kbox) 
  do while( iatom /= 0 )
    
    if ( skip(iatom) ) then
      iatom = inextbox(iatom)
      cycle
    end if

    d = (x-coor(iatom,1))**2 + &
        (y-coor(iatom,2))**2 + &
        (z-coor(iatom,3))**2 

    if ( d < vdwrad2 ) then
      d = dsqrt(d)
      dmin = dmin1(d,dmin)
      overlap = overlap + kvdw*( vdwrad - d )**2
    end if

    iatom = inextbox(iatom)
  end do

end subroutine overlapcell
