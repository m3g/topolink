!
! Determines if a atom is solvent-accessible
!
! An atom is considered solvent accessible here if there
! is a path of empty cells sharing faces (not vertices only)
! leading to the atom cell to the surface cell
!
subroutine solventaccess()

  use functionpars
  use linkedcells
  implicit none
  integer :: i, ibox, jbox, kbox
 
  ! 
  do i = 1, natom



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

