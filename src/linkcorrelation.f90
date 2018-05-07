!
! LinkCorrelation
!
! Computes the correlation of topological links given TopoLink log files.
!
! Leandro Martinez
! Institute of Chemistry - University of Campinas
! http://leandro.iqm.unicamp.br
!
! This package obtains the matrix of link correlations given
! a list of TopoLink output logs.
!
! For each link pair, the correlations output will be:
!
! type 0 : Print the fraction of structures that satisfy  
!          both links minus the fraction of structures that   
!          satisfy one link or the other.                 
!                                                        
! type 1 : Print the fraction of structures that satisfy
!          both links 
!
! type 2 : Print the fraction of structures that fail to
!          satisfy both links
!
! type 3 : Print the fraction of structures that satisfy
!          one OR the other link
!
! type 4 : Print the overall correlation, that is, the sum of the fraction
!          of structures that satify both links or fail to satisfy
!          both links, minus the fraction of structures that satisfy
!          one link OR the other.
!
! The default option is type = 1, which is easier to visualize.
!

program linkcorrelation

  use topolink_data
  use topolink_operations

  implicit none
  integer :: i, j, ilink, imodel, type
  integer :: nargs, nmodels, ioerr, maxlinks, nlinks
  character(len=200) :: loglist, record, line
  double precision, allocatable :: correlation(:,:), fraction(:)
  type(specific_link) :: linktemp
  type(modeldata), allocatable :: model(:)

  ! Read list of log files from the command line

  nargs = iargc()
  if ( nargs /= 1 .and. nargs /= 3 ) then
    write(*,*)
    write(*,*) ' Run with: linkcorrelation loglist.txt [-type type]'
    write(*,*)
    write(*,*) ' type 0 : Print the fraction of structures that satisfy '
    write(*,*) '          both links minus the fraction of structures that  '
    write(*,*) '          satisfy one link or the other.                '
    write(*,*) '                                                        '
    write(*,*) ' type 1 : Print the fraction of structures that satisfy '
    write(*,*) '          both links                                    '
    write(*,*) '                                                        '
    write(*,*) ' type 2 : Print the fraction of structures that fail to '
    write(*,*) '          satisfy both links                            '
    write(*,*) '                                                        '
    write(*,*) ' type 3 : Print the fraction of structures that satisfy '
    write(*,*) '          one OR the other link                         '
    write(*,*) '                                                        '
    write(*,*) ' type 4 : Print the overall correlation, that is, the sum of the fraction'
    write(*,*) '          of structures that satify both links or fail to satisfy        '
    write(*,*) '          both links, minus the fraction of structures that satisfy      '
    write(*,*) '          one link OR the other.                                         '
    write(*,*)
    stop
  end if
  call getarg(1,loglist)
  type = 1
  if ( nargs == 3 ) then
    call getarg(3,record)
    read(record,*,iostat=ioerr) type
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not read type parameter. '
      stop
    end if
    if ( type < 0 .or. type > 4 ) then
      write(*,*) ' ERROR: Invalid value for type parameter.'
      stop
    end if
  end if

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
    read(10,"(a200)",iostat=ioerr) record
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
    read(10,"(a200)",iostat=ioerr) record
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
      if ( line(4:11) == "RESULT0:") read(line(12:200),*,iostat=ioerr) model(imodel)%nobscons
      if ( ioerr /= 0 ) model(imodel)%nobscons = 0
      if ( line(4:11) == "RESULT1:") read(line(12:200),*,iostat=ioerr) model(imodel)%ntopcons
      if ( ioerr /= 0 ) model(imodel)%ntopcons = 0
      if ( line(4:11) == "RESULT2:") read(line(12:200),*,iostat=ioerr) model(imodel)%ntopnot
      if ( ioerr /= 0 ) model(imodel)%ntopnot = 0
      if ( line(4:11) == "RESULT3:") read(line(12:200),*,iostat=ioerr) model(imodel)%nmiss
      if ( ioerr /= 0 ) model(imodel)%nmiss = 0
      if ( line(4:11) == "RESULT4:") read(line(12:200),*,iostat=ioerr) model(imodel)%sumscores
      if ( ioerr /= 0 ) model(imodel)%sumscores = 0.
      if ( line(4:11) == "RESULT5:") read(line(12:200),*,iostat=ioerr) model(imodel)%likeli
      if ( ioerr /= 0 ) model(imodel)%likeli = 0.
      if ( line(4:11) == "RESULT6:") read(line(12:200),*,iostat=ioerr) model(imodel)%loglikeli
      if ( ioerr /= 0 ) model(imodel)%loglikeli = 0.
      if ( line(4:11) == "RESULT7:") read(line(12:200),*,iostat=ioerr) model(imodel)%usrlike
      if ( ioerr /= 0 ) model(imodel)%usrlike = 0.
      if ( line(4:11) == "RESULT8:") read(line(12:200),*,iostat=ioerr) model(imodel)%usrloglike
      if ( ioerr /= 0 ) model(imodel)%usrloglike = 0.
    end do
    close(20)
    model(imodel)%nlinks = ilink
  end do

  ! Indexing the links

  imodel = 1
  do i = 1, model(imodel)%nlinks
    model(imodel)%linkindex(i) = i
  end do
  do imodel = 2, nmodels 
    do i = 1, model(imodel)%nlinks
      do j = 1, model(1)%nlinks 
        if ( model(imodel)%link(i) .eq. model(1)%link(j) ) then
          model(imodel)%linkindex(i) = model(1)%linkindex(j)
          exit
        end if
      end do
    end do 
  end do

  ! Computing link satisfaction correlations

  write(*,"(a)") '# Computing link correlations: '
  write(*,"(a)") '#'
  write(*,"(a)") '# Reseting link correlation array ... '
  allocate(correlation(nlinks,nlinks),fraction(nlinks))
  do i = 1, nlinks
    fraction(i) = 0.d0
    do j = i, nlinks
      correlation(i,j) = 0.d0
    end do
  end do

  write(*,"(a)") '# Computing link correlations ... '
  do imodel = 1, nmodels
    do i = 1, model(imodel)%nlinks
      if ( ( model(imodel)%link(i)%status == 0 .or. &
             model(imodel)%link(i)%status == 1 .or. &
             model(imodel)%link(i)%status == 5 ) ) then
        fraction(i) = fraction(i) + 1.d0
      end if
      do j = i, model(imodel)%nlinks
        !
        ! Satisfied at the same time
        !
        if ( ( model(imodel)%link(i)%status == 0 .or. &
               model(imodel)%link(i)%status == 1 .or. &
               model(imodel)%link(i)%status == 5 ) .and. &
             ( model(imodel)%link(j)%status == 0 .or. &
               model(imodel)%link(j)%status == 1 .or. &
               model(imodel)%link(j)%status == 5 ) ) then
          if ( type == 0 .or. type == 1 .or. type == 4 ) then
            correlation(i,j) = correlation(i,j) + 1.d0
          end if
          cycle
        end if
        !
        ! Not satisfied at the same time
        !
        if ( ( model(imodel)%link(i)%status /= 0 .and. &
               model(imodel)%link(i)%status /= 1 .and. &
               model(imodel)%link(i)%status /= 5 ) .and. &
             ( model(imodel)%link(j)%status /= 0 .and. &
               model(imodel)%link(j)%status /= 1 .and. &
               model(imodel)%link(j)%status /= 5 ) ) then
          if ( type == 2 .or. type == 4 ) then 
            correlation(i,j) = correlation(i,j) + 1.d0
          end if
          cycle
        end if
        !
        ! One is satified, the other is not (if got here)
        !
        if ( type == 0 .or. type == 3 .or. type == 4 ) then
          correlation(i,j) = correlation(i,j) - 1.d0
        end if
      end do
    end do
  end do

  ! Writting the links

  do i = 1, nlinks
    fraction(i) = fraction(i) / dble(nmodels)
    write(*,"(a,i5,a,a,i5,a,f8.3)") &
          "#",i, trim(print_link(model(1)%link(i))), " (",i-1,") ", fraction(i) 
  end do

  do i = 1, nlinks
    do j = i, nlinks
      correlation(i,j) = correlation(i,j) / dble(nmodels)
      correlation(j,i) = correlation(i,j)
    end do
    write(*,"( 1000(tr1,f5.2) )") (correlation(i,j),j=1,nlinks)
  end do

end program linkcorrelation



