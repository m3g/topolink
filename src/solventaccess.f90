!
! subroutine solventacess: Determines if the residues are solvent-accessible
!
! An atom is considered solvent accessible here if there
! is a path of empty cells sharing faces (not vertices only)
! leading to the atom cell to the surface cell
!
! A residue is considered solvent-accessible if at least
! one of its atoms is solvent-accessible
!
! L. Martinez
! Institute of Chemistry - University of Campinas
! http://leandro.iqm.unicamp.br
! Nov 26, 2016
!
subroutine solventaccess(atom)

  use topolink_data
  use functionpars
  use linkedcells
  use ioformat
  use inputoptions, only : printaccessible
  implicit none
  integer :: ibox, jbox, kbox, n, i, j
  logical :: checkfaces
  logical, allocatable :: acc(:,:,:), aux(:,:,:)
  type(pdbatom) :: atom(natoms)

  ! Allocate logical accessible arrays

  allocate( acc(nboxesx,nboxesy,nboxesz), &
            aux(nboxesx,nboxesy,nboxesz) )

  ! Initialize all cells as unaccessible

  do ibox = 1, nboxesx
    do jbox = 1, nboxesy
      do kbox = 1, nboxesz
        acc(ibox,jbox,kbox) = .false.
      end do
    end do
  end do

  ! Cells at the borders of the grid are, of course, accessible

  ibox = 1
  do jbox = 1, nboxesy
    do kbox = 1, nboxesz
      acc(ibox,jbox,kbox) = .true.
    end do
  end do

  ibox = nboxesx
  do jbox = 1, nboxesy
    do kbox = 1, nboxesz
      acc(ibox,jbox,kbox) = .true.
    end do
  end do

  jbox = 1
  do ibox = 1, nboxesx
    do kbox = 1, nboxesz
      acc(ibox,jbox,kbox) = .true.
    end do
  end do

  jbox = nboxesy
  do ibox = 1, nboxesx
    do kbox = 1, nboxesz
      acc(ibox,jbox,kbox) = .true.
    end do
  end do

  kbox = 1
  do ibox = 1, nboxesx
    do jbox = 1, nboxesy
      acc(ibox,jbox,kbox) = .true.
    end do
  end do

  kbox = nboxesz
  do ibox = 1, nboxesx
    do jbox = 1, nboxesy
      acc(ibox,jbox,kbox) = .true.
    end do
  end do

  ! Empty cells sharing faces with accessible boxes are also accessible (iterate
  ! until convergence)

  n = 1
  do while( n > 0 ) 
    n = 0
    do ibox = 2, nboxesx - 1 
      do jbox = 2, nboxesy - 1
        do kbox = 2, nboxesz - 1
          if ( .not. acc(ibox,jbox,kbox ) ) then
            if ( ifirstbox(ibox,jbox,kbox) == 0 ) then
              if ( checkfaces(ibox,jbox,kbox,acc) ) then 
                acc(ibox,jbox,kbox) = .true.
                n = n + 1
              end if
            end if
          end if
        end do
      end do
    end do
  end do

  ! Finally, the cells containing atoms and sharing faces with 
  ! accessible cells are the important accessible cells

  do ibox = 2, nboxesx - 1 
    do jbox = 2, nboxesy - 1
      do kbox = 2, nboxesz - 1
        aux(ibox,jbox,kbox) = .false.
        if ( checkfaces(ibox,jbox,kbox,acc) ) aux(ibox,jbox,kbox) = .true.
      end do
    end do
  end do
  do ibox = 2, nboxesx - 1 
    do jbox = 2, nboxesy - 1
      do kbox = 2, nboxesz - 1
        if ( aux(ibox,jbox,kbox) ) acc(ibox,jbox,kbox) = .true.
      end do
    end do
  end do

  deallocate(aux)

  ! Check, residue by residue, which are solvent accessible

  do i = 1, natoms
    atom(i)%accessible = .false.
    atom(i)%residue%accessible = .false.
  end do
  n = 0
  do i = 1, natoms
    ibox = int( (coor(i,1)-xmin)/vdwrad ) + 1 
    jbox = int( (coor(i,2)-ymin)/vdwrad ) + 1 
    kbox = int( (coor(i,3)-zmin)/vdwrad ) + 1 
    if ( acc(ibox,jbox,kbox) ) then
      atom(i)%accessible = .true.
      if ( .not. atom(i)%residue%accessible ) then
        do j = atom(i)%residue%firstatom, atom(i)%residue%lastatom
          atom(j)%residue%accessible = .true.
        end do
      end if
      n = n + 1
    end if
  end do
  call writelog(blank)
  write(str,*) ' Number of atoms accessible to solvent: ', n ; call writelog(str)

  ! Count number of residues accessible to solvent

  n = 0
  i = 1
  do 
    i = atom(i)%residue%lastatom
    if ( atom(i)%residue%accessible ) n = n + 1
    if ( i == natoms ) exit  
    i = i + 1
  end do
  write(str,*) ' Number of residues accessible to solvent: ', n ; call writelog(str)

  ! Print residue accessibility to solvent

  if ( printaccessible ) then
    open(10,file="accessibility.pdb")
    do i = 1, natoms
      atom(i)%b = 1.
      if ( .not. atom(i)%residue%accessible ) atom(i)%b = 0.
      write(10,"(a)") trim(adjustl(print_pdbatom(atom(i))))
    end do
    ! Prints region protected from solvent
    !n = natoms
    !do ibox = 1, nboxesx
    !  do jbox = 1, nboxesy
    !    do kbox = 1, nboxesz
    !      if ( .not. acc(ibox,jbox,kbox) ) then
    !        n = n + 1
    !        writeatom%index = n
    !        writeatom%name = "O"
    !        writeatom%residue%name = "WAT"
    !        writeatom%residue%index = n
    !        writeatom%x = xmin + (ibox-1)*vdwrad + vdwrad/2.
    !        writeatom%y = ymin + (jbox-1)*vdwrad + vdwrad/2.
    !        writeatom%z = zmin + (kbox-1)*vdwrad + vdwrad/2.
    !        write(10,"(a)") trim(adjustl(print_pdbhetatm(writeatom)))
    !      end if
    !    end do
    !  end do
    !end do
    close(10)
  end if

  deallocate(acc)

end subroutine solventaccess

! 
! This function returns .true. if the box shares a face
! with an accessible box, and .false. otherwise
!

logical function checkfaces(ibox,jbox,kbox,acc)
  
  use linkedcells
  implicit none
  integer :: ibox, jbox, kbox
  logical :: acc(nboxesx,nboxesy,nboxesz)

  checkfaces = .false.
  if ( acc(ibox-1,jbox,kbox) .or. &
       acc(ibox,jbox-1,kbox) .or. &
       acc(ibox,jbox,kbox-1) .or. &
       acc(ibox+1,jbox,kbox) .or. &
       acc(ibox,jbox+1,kbox) .or. &
       acc(ibox,jbox,kbox+1) &
     ) then
    checkfaces = .true.
  end if

end function checkfaces
