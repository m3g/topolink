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
! type 4 : Print the the sum of the fraction
!          of structures that satify both links or fail to satisfy
!          both links, minus the fraction of structures that satisfy
!          one link OR the other.
!
! type 5 : Prints the Phi correlation, that is,
!          phi(1,2) = (n11*n00-n10*n01)/ sqrt(n1x*n0x*nx1*nx0) 
!          where: n11: number of models satisfying both links
!                 n00: number of models not satisfying any of the links
!                 n10: number of models satyisfying the first but not the second link
!                 n01: number of models satyisfying the second but not the first link
!                 n1x: number of models satyisfying the first link
!                 n0x: number of models not satyisfying the first link
!                 nx1: number of models satyisfying the second link
!                 nx0: number of models not satyisfying the second link
!
! The default option is type = 1, which is easier to visualize.
!

program linkcorrelation

  use ioformat, only : max_string_length, string_read
  use topolink_data
  use topolink_operations

  implicit none
  integer :: i, j, ilink, imodel, type
  integer :: nargs, nmodels, ioerr, maxlinks, nlinks
  character(len=max_string_length) :: loglist, record, line, format
  double precision, allocatable :: correlation(:,:), fraction(:)
  double precision, allocatable :: n11(:,:), n00(:,:), n10(:,:), n01(:,:), n1x(:,:), n0x(:,:), nx1(:,:), nx0(:,:)
  double precision :: phiden
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
    write(*,*) ' type 4 : Print the the sum of the fraction'
    write(*,*) '          of structures that satify both links or fail to satisfy        '
    write(*,*) '          both links, minus the fraction of structures that satisfy      '
    write(*,*) '          one link OR the other.                                         '
    write(*,*) '                                                        '
    write(*,*) ' type 5 : Prints the Phi correlation, that is,'
    write(*,*) '          phi(1,2) = (n11*n00-n10*n01)/ sqrt(n1x*n0x*nx1*nx0) '
    write(*,*) '          where: n11: number of models satisfying both links'
    write(*,*) '                 n00: number of models not satisfying any of the links'
    write(*,*) '                 n10: number of models satyisfying the first but not the second link'
    write(*,*) '                 n01: number of models satyisfying the second but not the first link'
    write(*,*) '                 n1x: number of models satyisfying the first link'
    write(*,*) '                 n0x: number of models not satyisfying the first link'
    write(*,*) '                 nx1: number of models satyisfying the second link'
    write(*,*) '                 nx0: number of models not satyisfying the second link'
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
    if ( type < 0 .or. type > 5 ) then
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
    read(10,string_read,iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    open(20,file=record,status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) cycle
    nmodels = nmodels + 1
    ! Check the number of links reported in this file
    nlinks = 0
    do 
      read(20,string_read,iostat=ioerr) line
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
    read(10,string_read,iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    open(20,file=record,status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) cycle
    imodel = imodel + 1
    model(imodel)%name = record
    ilink = 0
    do 
      read(20,string_read,iostat=ioerr) line
      if ( ioerr /= 0 ) exit
      if ( line(3:7) == "LINK:" ) then
        linktemp = read_link(line)
        if ( linktemp%observed ) then
          ilink = ilink + 1
          model(imodel)%link(ilink) = linktemp
        end if
      end if
      record = adjustl(line)
      if ( record(1:8) == "RESULT0:") read(record(9:max_string_length),*,iostat=ioerr) model(imodel)%nobscons
      if ( ioerr /= 0 ) model(imodel)%nobscons = 0
      if ( record(1:8) == "RESULT1:") read(record(9:max_string_length),*,iostat=ioerr) model(imodel)%ntopcons
      if ( ioerr /= 0 ) model(imodel)%ntopcons = 0
      if ( record(1:8) == "RESULT2:") read(record(9:max_string_length),*,iostat=ioerr) model(imodel)%ntopnot
      if ( ioerr /= 0 ) model(imodel)%ntopnot = 0
      if ( record(1:8) == "RESULT3:") read(record(9:max_string_length),*,iostat=ioerr) model(imodel)%nmiss
      if ( ioerr /= 0 ) model(imodel)%nmiss = 0
      if ( record(1:8) == "RESULT4:") read(record(9:max_string_length),*,iostat=ioerr) model(imodel)%sumscores
      if ( ioerr /= 0 ) model(imodel)%sumscores = 0.
      if ( record(1:8) == "RESULT5:") read(record(9:max_string_length),*,iostat=ioerr) model(imodel)%likeli
      if ( ioerr /= 0 ) model(imodel)%likeli = 0.
      if ( record(1:8) == "RESULT6:") read(record(9:max_string_length),*,iostat=ioerr) model(imodel)%loglikeli
      if ( ioerr /= 0 ) model(imodel)%loglikeli = 0.
      if ( record(1:8) == "RESULT7:") read(record(9:max_string_length),*,iostat=ioerr) model(imodel)%usrlike
      if ( ioerr /= 0 ) model(imodel)%usrlike = 0.
      if ( record(1:8) == "RESULT8:") read(record(9:max_string_length),*,iostat=ioerr) model(imodel)%usrloglike
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
  allocate(n11(nlinks,nlinks),&
           n00(nlinks,nlinks),&
           n10(nlinks,nlinks),&
           n01(nlinks,nlinks),&
           n1x(nlinks,nlinks),&
           nx1(nlinks,nlinks),&
           n0x(nlinks,nlinks),&
           nx0(nlinks,nlinks))
  do i = 1, nlinks
    fraction(i) = 0.d0
    do j = i, nlinks
      correlation(i,j) = 0.d0
      n11(i,j) = 0.
      n00(i,j) = 0.
      n10(i,j) = 0. 
      n01(i,j) = 0. 
      n1x(i,j) = 0. 
      nx1(i,j) = 0.
      n0x(i,j) = 0.
      nx0(i,j) = 0.
    end do
  end do

  write(*,"(a)") '# Computing link correlations ... '
  do imodel = 1, nmodels
    do i = 1, model(imodel)%nlinks
      do j = i, model(imodel)%nlinks
        !
        ! First satisfied
        !
        if ( model(imodel)%link(i)%status == 0 .or. &
             model(imodel)%link(i)%status == 1 .or. &
             model(imodel)%link(i)%status == 5 ) then
          ! Second satisfied
          if ( model(imodel)%link(j)%status == 0 .or. &
               model(imodel)%link(j)%status == 1 .or. &
               model(imodel)%link(j)%status == 5 ) then
            n11(i,j) = n11(i,j) + 1.
          ! Second not satisfied
          else
            n10(i,j) = n10(i,j) + 1
          end if
          cycle
        end if
        !
        ! First not satisfied
        !
        if ( model(imodel)%link(i)%status /= 0 .and. &
             model(imodel)%link(i)%status /= 1 .and. &
             model(imodel)%link(i)%status /= 5 ) then
          ! Second satisfied 
          if ( model(imodel)%link(j)%status == 0 .and. &
               model(imodel)%link(j)%status == 1 .and. &
               model(imodel)%link(j)%status == 5 ) then
            n01(i,j) = n01(i,j) + 1
          ! Second not satisfied
          else
            n00(i,j) = n00(i,j) + 1.
          end if
          cycle
        end if
      end do
    end do
  end do

  if ( type == 0 ) then
    do i = 1, nlinks
      do j = 1, nlinks
        correlation(i,j) = ( n11(i,j) - (n10(i,j) + n01(i,j)) ) / nmodels
      end do
    end do
  end if

  if ( type == 1 ) then
    do i = 1, nlinks
      do j = 1, nlinks
        correlation(i,j) = n11(i,j) / nmodels
      end do
    end do
  end if

  if ( type == 2 ) then
    do i = 1, nlinks
      do j = 1, nlinks
        correlation(i,j) = ( n10(i,j) + n01(i,j) + n00(i,j) ) / nmodels
        correlation(j,i) = correlation(i,j)
      end do
    end do
  end if

  if ( type == 3 ) then
    do i = 1, nlinks
      do j = 1, nlinks
        correlation(i,j) = ( n10(i,j) + n01(i,j) ) / nmodels
        correlation(j,i) = correlation(i,j)
      end do
    end do
  end if

  if ( type == 4 ) then
    do i = 1, nlinks
      do j = 1, nlinks
        correlation(i,j) = ( ( n11(i,j) + n00(i,j) ) - ( n10(i,j) + n01(i,j) ) ) / nmodels 
        correlation(j,i) = correlation(i,j)
      end do
    end do
  end if

  if ( type == 5 ) then
    do i = 1, nlinks
      do j = 1, nlinks
        n1x(i,j) = n11(i,j)+n10(i,j)
        n0x(i,j) = n01(i,j)+n00(i,j)
        nx1(i,j) = n11(i,j)+n01(i,j)
        nx0(i,j) = n10(i,j)+n00(i,j)
      end do
    end do
    do i = 1, nlinks
      do j = 1, nlinks
        phiden = n1x(i,j)*n0x(i,j)*nx1(i,j)*nx0(i,j) 
        if ( phiden > 0. ) then
          correlation(i,j) = ( n11(i,j)*n00(i,j) - n10(i,j)*n01(i,j) ) / sqrt( phiden ) 
          correlation(j,i) = correlation(i,j)
        else
          correlation(i,j) = 0.
          correlation(j,i) = correlation(i,j)
        end if
      end do
    end do
  end if

  ! Writting the links

  do i = 1, nlinks
    do j = 1, nlinks
      fraction(i) = n11(i,j) + n10(i,j) 
    end do
    fraction(i) = fraction(i) / nmodels
    write(*,"(a,i5,a,a,i5,a,f8.3)") &
          "#",i, trim(print_link(model(1)%link(i))), " (",i-1,") ", fraction(i) 
  end do

  write(format,*) "(",nlinks,"(tr1,f6.3))"
  do i = 1, nlinks
    write(*,format) (correlation(i,j),j=1,nlinks)
  end do

end program linkcorrelation



