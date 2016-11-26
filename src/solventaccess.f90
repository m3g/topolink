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
  logical :: checkfaces

integer :: j, n, nlast

  ! Allocate logical accessible array

  allocate( accessible(nboxesx,nboxesy,nboxesz) )

  ! Initialize all cells as unaccessible

  do ibox = 1, nboxesx
    do jbox = 1, nboxesy
      do kbox = 1, nboxesz
        accessible(ibox,jbox,kbox) = .false.
      end do
    end do
  end do

  ! Cells at the surface are solvent accessible

  ibox = 1
  do jbox = 1, nboxesy
    do kbox = 1, nboxesz
      accessible(ibox,jbox,kbox) = .true.
    end do
  end do

  ibox = nboxesx
  do jbox = 1, nboxesy
    do kbox = 1, nboxesz
      accessible(ibox,jbox,kbox) = .true.
    end do
  end do

  jbox = 1
  do ibox = 1, nboxesx
    do kbox = 1, nboxesz
      accessible(ibox,jbox,kbox) = .true.
    end do
  end do

  jbox = nboxesy
  do ibox = 1, nboxesx
    do kbox = 1, nboxesz
      accessible(ibox,jbox,kbox) = .true.
    end do
  end do

  kbox = 1
  do ibox = 1, nboxesx
    do jbox = 1, nboxesy
      accessible(ibox,jbox,kbox) = .true.
    end do
  end do

  kbox = nboxesz
  do ibox = 1, nboxesx
    do jbox = 1, nboxesy
      accessible(ibox,jbox,kbox) = .true.
    end do
  end do

  ! Empty cells vicinal to accessible boxes are accessible

  n = 0
  nlast = -1
  do while( n /= nlast ) 
    nlast = n
    n = 0
    do ibox = 2, nboxesx - 1 
      do jbox = 2, nboxesy - 1
        do kbox = 2, nboxesz - 1
          if ( ifirstbox(ibox,jbox,kbox) == 0 ) then
            accessible(ibox,jbox,kbox) = checkfaces(ibox,jbox,kbox)
            if( .not. accessible(ibox,jbox,kbox) ) then
              n = n + 1
            end if
          end if
        end do
      end do
    end do
write(*,*) n
  end do

  ! Is there any empty box not solvent-accessible remaining?

  i = 0
  do ibox = 1, nboxesx
    do jbox = 1, nboxesy
      do kbox = 1, nboxesz
        if( ifirstbox(ibox,jbox,kbox) == 0 .and. &
            .not. accessible(ibox,jbox,kbox) ) then
          i = i + 1
          write(*,*) 'O', (xmin+ibox*vdwrad+vdwrad/2), &
                          (ymin+jbox*vdwrad+vdwrad/2), &
                          (zmin+kbox*vdwrad+vdwrad/2)
        end if
      end do
    end do  
  end do
stop


end subroutine solventaccess

!
! This subroutine runs over cells in a wall and checks
! if this cell is empty, if it is, and it shares a face
! with an accessible cell, it is also accessible 
!

subroutine checkwall(axis,i)
              
  use linkedcells
  implicit none
  integer :: i, ibox, jbox, kbox
  character :: axis
  logical :: checkfaces

  if ( axis == "x" ) then
    ibox = i
    do jbox = i, nboxesy - i + 1
      do kbox = i, nboxesz - i + 1
        if ( ifirstbox(ibox,jbox,kbox) == 0 ) then
          accessible(ibox,jbox,kbox) = checkfaces(ibox,jbox,kbox)
        end if
      end do
    end do
    ibox = nboxesx - i + 1
    do jbox = i, nboxesy - i + 1
      do kbox = i, nboxesz - i + 1
        if ( ifirstbox(ibox,jbox,kbox) == 0 ) then
          accessible(ibox,jbox,kbox) = checkfaces(ibox,jbox,kbox)
        end if
      end do
    end do
  end if

  if ( axis == "y" ) then
    jbox = i
    do ibox = i, nboxesx - i + 1
      do kbox = i, nboxesz - i + 1
        if ( ifirstbox(ibox,jbox,kbox) == 0 ) then
          accessible(ibox,jbox,kbox) = checkfaces(ibox,jbox,kbox)
        end if
      end do
    end do
    jbox = nboxesy - i + 1
    do ibox = i, nboxesx - i + 1
      do kbox = i, nboxesz - i + 1
        if ( ifirstbox(ibox,jbox,kbox) == 0 ) then
          accessible(ibox,jbox,kbox) = checkfaces(ibox,jbox,kbox)
        end if
      end do
    end do
  end if

  if ( axis == "z" ) then
    kbox = i
    do ibox = i, nboxesx - i + 1
      do jbox = i, nboxesy - i + 1
        if ( ifirstbox(ibox,jbox,kbox) == 0 ) then
          accessible(ibox,jbox,kbox) = checkfaces(ibox,jbox,kbox)
        end if
      end do
    end do
    kbox = nboxesz - i + 1
    do ibox = i, nboxesx - i + 1
      do jbox = i, nboxesy - i + 1
        if ( ifirstbox(ibox,jbox,kbox) == 0 ) then
          accessible(ibox,jbox,kbox) = checkfaces(ibox,jbox,kbox)
        end if
      end do
    end do
  end if

end subroutine checkwall

! 
! This function returns .true. if the box shares a face
! with an accessible box, and .false. otherwise
!

logical function checkfaces(ibox,jbox,kbox)
  
  use linkedcells
  implicit none
  integer :: ibox, jbox, kbox

  if ( accessible(ibox-1,jbox,kbox) .or. &
    accessible(ibox,jbox-1,kbox) .or. &
    accessible(ibox,jbox,kbox-1) .or. &
    accessible(ibox+1,jbox,kbox) .or. &
    accessible(ibox,jbox+1,kbox) .or. &
    accessible(ibox,jbox,kbox+1) ) then
    checkfaces = .true.
  end if

end function checkfaces

