
program filtermodels

  use topolink_data
  use topolink_operations
  use flashsort

  implicit none

  !
  ! Model data type
  !

  type modeldata
    ! Model log file name
    character(len=200) :: name
    ! Number of links in this file
    integer :: nlinks
    ! Link results for this model
    type(specific_link), allocatable :: link(:)
    ! Results
    integer :: nobscons  ! RESULT0: Number of observations that are consistent with the structure.
    integer :: ntopcons  ! RESULT1: Number of topological distances consistent with all observations.
    integer :: ntopnot   ! RESULT2: Number of topological distances NOT consistent with observations.
    integer :: nmiss     ! RESULT3: Number of missing links.
    double precision :: sumscores  ! RESULT4: Sum of scores of observed links of all experiments.
    double precision :: likely     ! RESULT5: Likelyhood of the set of experimental results.
    double precision :: loglikely  ! RESULT6: Log-likelyhood of the set of experimental results.
    double precision :: usrlike    ! RESULT7: Likelyhood of the set of experimental results. (with user pgood and pbad)
    double precision :: usrloglike ! RESULT8: Log-likelyhood of the set of experimental results. (with user pgood and pbad)
    ! Index of link in overall model links lists
    integer, allocatable :: linkindex(:)
  end type modeldata

  integer :: i, j, ilink, imodel
  integer :: nargs, nmodels, ioerr, maxlinks, nlinks, ncut, ngood
  character(len=200) :: loglist, linklist, record, line
  type(specific_link) :: linktemp
  type(observed_link), allocatable :: link(:)
  type(modeldata), allocatable :: model(:)

  ! Read list of log files from the command line

  nargs = iargc()
  if ( nargs /= 3 ) then
    write(*,*) ' Run with: filtermodels loglist.txt linklist.txt [ncut] '
    stop
  end if
  call getarg(1,loglist)
  call getarg(2,linklist)
  call getarg(3,record)
  read(record,*) ncut

  ! Checks how many log files are available

  write(*,"(a)") '#'
  write(*,"(a)") '# Reading input file ... '
  write(*,"(a)") '#'
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
      if ( line(3:7) == "LINK:" ) then
        linktemp = read_link(line)
        if ( linktemp%observed ) then
          nlinks = nlinks + 1
        end if
      end if
    end do
    maxlinks = max0(nlinks,maxlinks)
    close(20)
  end do
  close(20)
  write(*,"(a,i5)") '# Number of topolink log files found: ', nmodels
  write(*,"(a,i5)") '# Maximum number of links in a file: ', maxlinks
  rewind(10)

  ! Read model data

  write(*,"(a)") '#'
  write(*,"(a)") '# Reading model data file ... '
  write(*,"(a)") '#'
  allocate(model(nmodels))
  do imodel = 1, nmodels
    allocate(model(imodel)%link(maxlinks),model(imodel)%linkindex(maxlinks))
  end do
  imodel = 0
  do
    read(10,*,iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    open(20,file=record,status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) cycle
    imodel = imodel + 1
    model(imodel)%name = record
    ilink = 0
    do 
      read(20,"(a200)",iostat=ioerr) line
      if ( ioerr /= 0 ) exit
      if ( line(3:7) == "LINK:" ) then
        linktemp = read_link(line)
        if ( linktemp%observed ) then
          ilink = ilink + 1
          model(imodel)%link(ilink) = linktemp
        end if
      end if
      if ( line(4:11) == "RESULT0:") read(line(12:200),*) model(imodel)%nobscons
      if ( line(4:11) == "RESULT1:") read(line(12:200),*) model(imodel)%ntopcons
      if ( line(4:11) == "RESULT2:") read(line(12:200),*) model(imodel)%ntopnot
      if ( line(4:11) == "RESULT3:") read(line(12:200),*) model(imodel)%nmiss
      if ( line(4:11) == "RESULT4:") read(line(12:200),*) model(imodel)%sumscores
      if ( line(4:11) == "RESULT5:") read(line(12:200),*) model(imodel)%likely
      if ( line(4:11) == "RESULT6:") read(line(12:200),*) model(imodel)%loglikely
      if ( line(4:11) == "RESULT7:") read(line(12:200),*) model(imodel)%usrlike
      if ( line(4:11) == "RESULT8:") read(line(12:200),*) model(imodel)%usrloglike
    end do
    close(20)
    model(imodel)%nlinks = ilink
  end do
  close(10)

  ! Writting the links

  do i = 1, nlinks
    write(*,"(a,i5,a,a,i5,a)") "#",i, trim(print_link(model(1)%link(i))), " (",i-1,")"
  end do

  ! Reading link list

  open(10,file=linklist,status='old',action='read',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open link list: ', trim(adjustl(linklist))
    stop
  end if
  nlinks = 0
  do
    read(10,"(a200)",iostat=ioerr) line
    if( ioerr /= 0 ) exit
    nlinks = nlinks + 1
  end do
  if ( nlinks == 0 ) then
    write(*,*) ' ERROR: Could not find any link in linklist. '
    stop
  end if
  rewind(10)
  allocate(link(nlinks))
  write(*,"(a,i5)") '# Number of links to be considered: ', nlinks
  nlinks = 0
  write(*,"(a)") '# Links to be considered: '
  do
    read(10,"(a200)",iostat=ioerr) line
    if( ioerr /= 0 ) exit
    nlinks = nlinks + 1
    read(line,*,iostat=ioerr) link(nlinks)%residue1%name, &
                              link(nlinks)%residue1%chain, &
                              link(nlinks)%residue1%index, &
                              link(nlinks)%residue2%name, &
                              link(nlinks)%residue2%chain, &
                              link(nlinks)%residue2%index
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not read link in line: ', trim(adjustl(line))
    end if
    write(*,"(a,i5,2(tr1,a),i5,2(tr1,a),i5)") '#', nlinks, &
                          link(nlinks)%residue1%name, &
                          link(nlinks)%residue1%chain, &
                          link(nlinks)%residue1%index, &
                          link(nlinks)%residue2%name, &
                          link(nlinks)%residue2%chain, &
                          link(nlinks)%residue2%index
  end do
  close(10)

  ! Finding which models satisfy all of the selected links
  
  do imodel = 1, nmodels
    ngood = 0
    linksdo : do i = 1, nlinks
      do j = 1, model(imodel)%nlinks
        if ( model(imodel)%link(j) .matches. link(i) ) then
          if ( model(imodel)%link(j)%status == 0 .or. &
               model(imodel)%link(j)%status == 1 .or. &
               model(imodel)%link(j)%status == 5 ) then
            ngood = ngood + 1
          end if
          cycle linksdo
        end if
      end do
    end do linksdo
    if ( ngood >= ncut ) then
      write(*,*) trim(adjustl(model(imodel)%name)), ngood
    end if
  end do

end program filtermodels


