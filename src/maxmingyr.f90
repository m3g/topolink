
module datatypes

  implicit none

  type pdbatom
    character(len=4) :: type
    double precision :: x, y, z
    logical :: error
  end type pdbatom
  
  type pdbmodel
     character(len=200) :: file
     integer :: natoms
     type(pdbatom), allocatable :: atom(:)
     double precision :: gyration
  end type pdbmodel

end module datatypes

program maxmingyr

  use datatypes
  implicit none
  integer :: ioerr, maxnatoms, i
  character(len=200) :: filelist, line
  type(pdbatom) :: atom, readatom
  type(pdbmodel) :: model, mingyrmodel, maxgyrmodel
  double precision :: xavg, yavg, zavg

  ! Read file list from command line

  call getarg(1,filelist)

  open(10,file=filelist,status='old',action='read',iostat=ioerr)
  maxnatoms = 0
  mingyrmodel%gyration = 1.d30
  maxgyrmodel%gyration = -1.d0
  do
    read(10,"(a200)",iostat=ioerr) model%file
    if( ioerr /= 0 ) exit
    open(20,file=model%file,status='old',action='read',iostat=ioerr)
    
    ! Read the number of atoms in this file

    model%natoms = 0
    do
      read(20,"(a200)",iostat=ioerr) line
      if ( ioerr /= 0 ) exit
      atom = readatom(line)
      if ( .not. atom%error .and. atom%type == "CA" ) model%natoms = model%natoms + 1
    end do 

    ! Allocate vector containing coordinates

    if ( model%natoms > maxnatoms ) then
      if ( allocated(model%atom) ) deallocate(model%atom)
      allocate(model%atom(model%natoms))
    end if

    ! Read atomic coordinates

    rewind(20)
    i = 0
    do
      read(20,"(a200)",iostat=ioerr) line
      if ( ioerr /= 0 ) exit
      atom = readatom(line)
      if ( .not. atom%error .and. atom%type == "CA" ) then
        i = i + 1
        model%atom(i) = atom
      end if
    end do
    close(20)

    ! Compute average coordinates
      
    xavg = 0.d0
    yavg = 0.d0
    zavg = 0.d0
    do i = 1, model%natoms 
      xavg = xavg + model%atom(i)%x
      yavg = yavg + model%atom(i)%y
      zavg = zavg + model%atom(i)%z
    end do
    xavg = xavg / model%natoms
    yavg = yavg / model%natoms
    zavg = zavg / model%natoms

    ! Compute gyration radius

    model%gyration = 0.d0
    do i = 1, model%natoms
      model%gyration = model%gyration + &
                 ( model%atom(i)%x - xavg )**2 + &
                 ( model%atom(i)%y - yavg )**2 + &
                 ( model%atom(i)%z - zavg )**2 
    end do
    model%gyration = dsqrt(model%gyration/model%natoms)

    if ( model%gyration < mingyrmodel%gyration ) then
      mingyrmodel%gyration = model%gyration
      mingyrmodel%file = model%file
    end if
    if ( model%gyration > maxgyrmodel%gyration ) then
      maxgyrmodel%gyration = model%gyration
      maxgyrmodel%file = model%file
    end if

  end do
  close(10)
  write(*,*) trim(adjustl(mingyrmodel%file)), mingyrmodel%gyration
  write(*,*) trim(adjustl(maxgyrmodel%file)), maxgyrmodel%gyration
  
end program maxmingyr

function readatom(line)

  use datatypes, only : pdbatom 
  implicit none
  integer :: ioerr
  character(len=200) :: line
  type(pdbatom) :: readatom

  readatom%error = .false.
  if ( line(1:4) == "ATOM" ) then
    read(line(13:16),*,iostat=ioerr) readatom%type
    if ( ioerr /= 0 ) then
      readatom%error = .true.
      return
    end if
    read(line(31:38),*,iostat=ioerr) readatom%x
    if( ioerr /= 0 ) then
      readatom%error = .true.
      return
    end if
    read(line(39:46),*,iostat=ioerr) readatom%y
    if( ioerr /= 0 ) then
      readatom%error = .true.
      return
    end if
    read(line(47:54),*,iostat=ioerr) readatom%z
    if( ioerr /= 0 ) then
      readatom%error = .true.
      return
    end if
  end if

end function readatom



