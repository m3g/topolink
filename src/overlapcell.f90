!
! Subroutine that computes the overlap in a cell
!

subroutine overlapcell(i,x,y,z,ibox,jbox,kbox,overlap)

  use functionpars
  use linkedcells
  implicit none
  integer :: ibox, jbox, kbox, iatom, i
  double precision :: x, y, z, overlap, d, vdwrad2

  if ( ibox < 1 .or. jbox < 1 .or. kbox < 1 ) return
  if ( ibox > nboxesx .or. jbox > nboxesy .or. kbox > nboxesz ) return
  
  vdwrad2 = sigma(i)**2
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
      dmin_maxviol = min(dmin_maxviol,sigma(i)-d)
      overlap = overlap + kvdw*( sigma(i) - d )**2
    end if

    iatom = inextbox(iatom)
  end do

end subroutine overlapcell
