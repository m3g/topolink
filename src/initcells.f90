!
! Subroutine that initializes the linkes cells 
!

subroutine initcells()

  use functionpars
  use linkedcells
  implicit none
  integer :: i, ibox, jbox, kbox

  ! Compute maximum and minimum coordinates of the structure

  xmin = coor(1,1)
  ymin = coor(1,2)
  zmin = coor(1,3)
  xmax = coor(1,1)
  ymax = coor(1,2)
  zmax = coor(1,3)
  do i = 2, natoms
    xmin = min(xmin,coor(i,1))
    ymin = min(ymin,coor(i,2))
    zmin = min(zmin,coor(i,3))
    xmax = max(xmax,coor(i,1))
    ymax = max(ymax,coor(i,2))
    zmax = max(zmax,coor(i,3))
  end do

  ! The maximum and minimum reference coordinates have added vdwrad distance

  xmin = xmin - vdwrad
  ymin = ymin - vdwrad
  zmin = zmin - vdwrad
  xmax = xmax + vdwrad
  ymax = ymax + vdwrad
  zmax = zmax + vdwrad

  ! Compute the number of boxes

  nboxesx = int( ( xmax - xmin ) / vdwrad ) + 1
  nboxesy = int( ( ymax - ymin ) / vdwrad ) + 1
  nboxesz = int( ( zmax - zmin ) / vdwrad ) + 1

  ! Allocate necessary arrays

  allocate( ifirstbox(nboxesx,nboxesy,nboxesz), inextbox(natoms) )
  do ibox = 1, nboxesx
    do jbox = 1, nboxesy
      do kbox = 1, nboxesz
        ifirstbox(ibox,jbox,kbox) = 0
      end do
    end do
  end do

  ! Add structure atoms to its boxes (permanently, as they don't move here)

  do i = 1, natoms
    ibox = int( (coor(i,1)-xmin)/vdwrad ) + 1 
    jbox = int( (coor(i,2)-ymin)/vdwrad ) + 1 
    kbox = int( (coor(i,3)-zmin)/vdwrad ) + 1 
    inextbox(i) = ifirstbox(ibox,jbox,kbox)
    ifirstbox(ibox,jbox,kbox) = i
  end do

end subroutine initcells

