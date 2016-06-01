
program evalmodels

  use topolink_data
  use topolink_operations
  use flashsort

  implicit none

  type resultlist
    integer :: r0, r1, r2, r3
    double precision :: r4, r5, r6, r7, r8
  end type resultlist

  integer :: i
  integer :: nargs, nmodels, ioerr, imodel, ilink, maxlinks, nlinks
  double precision :: scoremin
  character(len=200) :: loglist, record, record2, line
  character(len=200), allocatable :: name(:)
 
  type(resultlist), allocatable :: result(:)
  double precision, allocatable :: score(:), order(:)

  type(specific_link) :: linktemp
  type(specific_link), allocatable :: link(:,:)

  ! Read list of log files from the command line

  nargs = iargc()
  if ( nargs /= 1 ) then
    write(*,*) ' Run with: evalmodels loglist.txt '
    stop
  end if
  call getarg(1,loglist)

  ! Checks how many log files are available

  write(*,*) '#'
  write(*,*) '# Reading input file ... '
  write(*,*) '#'
  open(10,file=loglist,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open log list file. '
    stop
  end if
  nmodels = 0
  maxlinks = 0
  do 
    read(10,*,iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    open(20,file=record,status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) cycle
    nmodels = nmodels + 1
    ! Check the number of links reported in this file
    nlinks = 0
    do 
      read(20,"(a200)",iostat=ioerr) line
      if ( ioerr /= 0 ) exit
      if ( line(3:7) == "LINK:" ) nlinks = nlinks + 1
    end do
    maxlinks = max0(nlinks,maxlinks)
    close(20)
  end do
  close(20)
  write(*,*) '# Number of topolink log files found: ', nmodels
  write(*,*) '# Maximum number of links in a file: ', maxlinks
  rewind(10)

  ! Read model data

  write(*,*) '#'
  write(*,*) '# Reading model data file ... '
  write(*,*) '#'
  allocate(score(nmodels),result(nmodels))
  allocate(name(nmodels),link(nmodels,maxlinks))
  imodel = 0
  do
    read(10,*,iostat=ioerr) record, record2
    if ( ioerr /= 0 ) exit
    open(20,file=record,status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) cycle
    imodel = imodel + 1
    name(imodel) = record
    ilink = 0
    do 
      read(20,"(a200)",iostat=ioerr) line
      if ( ioerr /= 0 ) exit
      if ( line(3:7) == "LINK:" ) then
        linktemp = read_link(line)
        ilink = ilink + 1
        link(imodel,ilink) = linktemp
      end if
      if ( line(4:11) == "RESULT0:") read(line(12:200),*) result(imodel)%r0
      if ( line(4:11) == "RESULT1:") read(line(12:200),*) result(imodel)%r1
      if ( line(4:11) == "RESULT2:") read(line(12:200),*) result(imodel)%r2
      if ( line(4:11) == "RESULT3:") read(line(12:200),*) result(imodel)%r3
      if ( line(4:11) == "RESULT4:") read(line(12:200),*) result(imodel)%r4
      if ( line(4:11) == "RESULT5:") read(line(12:200),*) result(imodel)%r5
      if ( line(4:11) == "RESULT6:") read(line(12:200),*) result(imodel)%r6
      if ( line(4:11) == "RESULT7:") read(line(12:200),*) result(imodel)%r7
      if ( line(4:11) == "RESULT8:") read(line(12:200),*) result(imodel)%r8
    end do
    close(20)
    read(record2,*,iostat=ioerr) score(imodel)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not read score of model: ', trim(adjustl(record))
    end if
  end do

  ! Ordering the models in terms of score

  scoremin = 0.d0
  do imodel = 1, nmodels
    scoremin = dmin1(scoremin,score(imodel))
  end do

  mflash = 1 + nmodels/10
  allocate(indflash(nmodels),lflash(nmodels),order(nmodels))
  do imodel = 1, nmodels
    !order(imodel) = score(imodel)
    order(imodel) = result(imodel)%r5*(exp(-1.d0*(score(imodel)-scoremin)/0.593d0))
  end do
  call flash1(order,nmodels,lflash,mflash,indflash)
 
  ! Write something

  do imodel = 1, nmodels
    i = indflash(imodel)
    write(*,*) score(i), result(i)%r5, order(imodel), trim(adjustl(name(i)))
  end do

    
  



  


end program evalmodels



