!
! Program CompareCols
!
! This is a simple program that reads scores from two different
! files containing model names to specify to which model the score
! refers, and writes one score as a function of the other, ordering
! the models by one of the scores
!
! L. Martinez
! Instiute of Chemistry - University of Campinas
! September 13, 2016
! http://leandro.iqm.unicamp.br/topolink
!

module comparecols_data

  use ioformat, only : max_string_length
  type data
    character(len=max_string_length) :: name
    double precision :: score(2)
    logical :: foundpair
  end type data

end module comparecols_data

program comparecols

  use ioformat 
  use string_operations
  use comparecols_data
  implicit none
  integer :: imodel, i, ioerr, npairs
  integer :: nargs, nmodels, scol1, scol2, ncol1, ncol2, sort
  double precision :: score, scoremax(2), scoremin(2)
  character(len=max_string_length) :: file1, file2, output, record, string, name
  integer :: model_index
  logical :: error
  type(data), allocatable :: model(:)

  ! Print title

  call title()
  write(*,*) ' COMPARECOLS: Compares two scores in general score files specified by model names. '
  write(*,*)
  write(*,dashes)

  ! Read list of log files from the command line

  nargs = iargc()
  if ( nargs /= 8 ) then
    write(*,*)
    write(*,*) ' ERROR: Run with: '
    write(*,*) ' comparecols score1.dat [score col] [name col] scores2.dat [score col] [name col] -s[1/2] output.dat'
    write(*,*)
    write(*,*) ' Where: score1.dat and score2.dat are the data files containing the scores.'
    write(*,*) 
    write(*,*) ' [score col] is an integer number indicating the column of the score in each file. ' 
    write(*,*) ' [name col] is an integer number indicating the column of the model name in each file. ' 
    write(*,*) ' -s0: The output will be sorted by model name. ' 
    write(*,*) ' -s1: The output will be sorted according to the score of the first file, ascending order. ' 
    write(*,*) ' -s2: The output will be sorted according to the score of the second file, ascending order. ' 
    write(*,*) ' -s3: The output will be sorted according to the score of the first file, descending order. ' 
    write(*,*) ' -s4: The output will be sorted according to the score of the second file, descending order. ' 
    write(*,*)
    write(*,*) ' output.dat is the output file. '
    write(*,*)
    write(*,*) ' More details at: http://leandro.iqm.unicamp/topolink '
    write(*,*)
    write(*,hashes)
    stop
  end if

  !
  ! First file
  !

  call getarg(1,file1)

  call getarg(2,record)
  read(record,*,iostat=ioerr) scol1
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read score column of first file (second argument). '
    stop
  end if
  call getarg(3,record)
  read(record,*,iostat=ioerr) ncol1
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read name column of first file (third argument). '
    stop
  end if

  !
  ! Second file
  !

  call getarg(4,file2)

  call getarg(5,record)
  read(record,*,iostat=ioerr) scol2
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read score column of second file (fifth argument). '
    stop
  end if
  call getarg(6,record)
  read(record,*,iostat=ioerr) ncol2
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read score column of second file (sixth argument). '
    stop
  end if

  ! Read from which column the scores have to be read from the scorelist

  sort = 0
  call getarg(7,record)
  read(record(3:length(record)),*,iostat=ioerr) sort
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read sort parameter -s0, -s1, -s2, -s3, -s4'
    stop
  end if

  ! Output file

  call getarg(8,output)

  ! Open first score file

  open(10,file=file1,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(file1))
    stop
  end if
  nmodels = 0
  do
    read(10,string_read,iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    read(record,*,iostat=ioerr) (string, i = 1, scol1)
    if ( ioerr /= 0 ) cycle
    read(string,*,iostat=ioerr) score
    if ( ioerr /= 0 ) cycle
    read(record,*,iostat=ioerr) (name, i = 1, ncol1)
    if ( ioerr /= 0 ) cycle
    nmodels = nmodels + 1
  end do
  write(*,*) ' Number of models found in the first file: ', nmodels

  allocate(model(nmodels))

  ! Reading data from file 1
  
  rewind(10)
  imodel = 0
  do
    read(10,string_read,iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    backspace(10) ; read(10,*,iostat=ioerr) (string, i = 1, scol1)
    if ( ioerr /= 0 ) cycle
    read(string,*,iostat=ioerr) score
    if ( ioerr /= 0 ) cycle
    backspace(10) ; read(10,*,iostat=ioerr) (name, i = 1, ncol1)
    if ( ioerr /= 0 ) cycle
    imodel = imodel + 1
    model(imodel)%name = basename(name)
    model(imodel)%score(1) = score
  end do
  close(10)

  ! Sorting models by name

  call sort_by_name(nmodels,model)

  ! Reading data from file 2

  open(10,file=file2,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(file2))
    stop
  end if
  do imodel = 1, nmodels
    model(imodel)%foundpair = .false.
  end do
  npairs = 0
  do
    read(10,string_read,iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    backspace(10) ; read(10,*,iostat=ioerr) (string, i = 1, scol2)
    if ( ioerr /= 0 ) cycle
    read(string,*,iostat=ioerr) score
    if ( ioerr /= 0 ) cycle
    backspace(10) ; read(10,*,iostat=ioerr) (name, i = 1, ncol2)
    if ( ioerr /= 0 ) cycle
    name = basename(name)
    imodel = model_index(name,model,nmodels,error)
    if ( .not. error ) then
      model(imodel)%foundpair = .true.
      model(imodel)%score(2) = score
      npairs = npairs + 1
    end if
  end do
  close(10)
  write(*,*) ' Number of pairs found: ', npairs
  if ( npairs == 0 ) then
    write(*,*) ' ERROR: No corresponding names were found between the the two files. '
    write(*,*) '        Probably the name column set for one of the files is not correct. '
    stop
  end if

  ! Sorting models according to the user choice

  if ( sort == 1 .or. sort == 2 .or. sort == 3 .or. sort == 4 ) then
    call sort_by_score(nmodels,model,sort)
  end if

  open(10,file=output)
  write(10,"(a,a)") "# File 1: ", trim(file1)
  write(10,"(a,i5)") "# Score column 1: ", scol1
  write(10,"(a,i5)") "# Name column 1: ", ncol1
  write(10,"(a,a)") "# File 2: ", trim(file2)
  write(10,"(a,i5)") "# Score column 2: ", scol2
  write(10,"(a,i5)") "# Name column 2: ", ncol2 
  write(10,"(a)") "#"
  if ( sort == 0 ) then
    write(10,"(a)") "# Models sorted by name."
  else
    if ( sort == 1 ) write(10,"(a,i2)") "# Models sorted by score 1 in ascending order. "
    if ( sort == 2 ) write(10,"(a,i2)") "# Models sorted by score 2 in ascending order. "
    if ( sort == 3 ) write(10,"(a,i2)") "# Models sorted by score 1 in descending order. "
    if ( sort == 4 ) write(10,"(a,i2)") "# Models sorted by score 2 in descending order. "
  end if
  write(10,"(a)") "#"
  write(10,"(a)") "# SCORE1: Score read from file 1."
  write(10,"(a)") "# SCORE2: Score read from file 2."
  write(10,"(a)") "#"
  write(10,"(a)") "# MIN1: Minimum score 1 found up to this model."
  write(10,"(a)") "# MAX1: Maximum score 1 found up to this model."
  write(10,"(a)") "# MIN2: Minimum score 2 found up to this model."
  write(10,"(a)") "# MAX2: Maximum score 2 found up to this model."
  write(10,"(a)") "#"
  write(10,"('#     SCORE1        SCORE2          MIN1          MAX1          MIN2          MAX2        NAME')")
  scoremax(1) = -1.d30
  scoremin(1) =  1.d30
  scoremax(2) = -1.d30
  scoremin(2) =  1.d30
  do imodel = 1, nmodels
    if ( model(imodel)%foundpair ) then
      scoremax(1) = max(scoremax(1),model(imodel)%score(1))
      scoremin(1) = min(scoremin(1),model(imodel)%score(1))
      scoremax(2) = max(scoremax(2),model(imodel)%score(2))
      scoremin(2) = min(scoremin(2),model(imodel)%score(2))
      write(10,"( f12.5, 5(tr2, f12.5), tr2, a )") &
               model(imodel)%score(1), model(imodel)%score(2),&
               scoremin(1), scoremax(1), scoremin(2), scoremax(2),&      
               trim(adjustl(model(imodel)%name))
    end if
  end do
  close(10)
  write(*,*)
  write(*,*) ' Wrote output file: ', trim(adjustl(output))
  write(*,*)
  write(*,hashes)

end program comparecols

!
! Function that determines the index of a model from its name
!

function model_index(name,model,n,error)
 
  use ioformat, only : max_string_length
  use comparecols_data
  implicit none
  integer :: model_index
  integer :: n, imax, imin, iavg
  character(len=max_string_length) :: name
  logical :: error
  type(data) :: model(n)

  imin = 1
  imax = n
  error = .false.
  if ( name < model(1)%name ) error = .true.
  if ( name > model(n)%name ) error = .true.
  if ( .not. error ) then
    do
      if ( name == model(imin)%name ) then 
        model_index = imin
        return
      end if
      if ( name == model(imax)%name ) then
        model_index = imax
        return
      end if
      if ( imax == imin ) then
        error = .true.
        exit
      end if
      iavg = imin + ( imax - imin ) / 2
      if ( name >= model(iavg)%name ) imin = iavg + 1
      if ( name <= model(iavg)%name ) imax = iavg
    end do
  end if
  if ( error ) then
    write(*,*) ' WARNING: A model is listed in a fist but not second file: ', trim(adjustl(name))
  end if

end function model_index

!
! Subroutines to sort models
!

!
! Sort by names
!

subroutine sort_by_name(n,model)

  use comparecols_data
  implicit none
  integer :: i, j, n
  type(data) :: model(n), modeltemp

  do i = 1, n-1
    j = i + 1
    do while( model(j-1)%name > model(j)%name )
      modeltemp = model(j-1)
      model(j-1) = model(j) 
      model(j) = modeltemp
      j = j - 1
      if ( j == 1 ) exit
    end do
  end do

end subroutine sort_by_name

!
! Sort by score
!

subroutine sort_by_score(n,model,sort)

  use comparecols_data
  implicit none
  integer :: i, j, n, iscore, sort
  type(data) :: model(n), modeltemp

  if ( sort == 1 .or. sort == 2 ) then
    iscore = sort
    do i = 1, n-1
      j = i + 1
      do while( model(j-1)%score(iscore) > model(j)%score(iscore) )
        modeltemp = model(j-1)
        model(j-1) = model(j) 
        model(j) = modeltemp
        j = j - 1
        if ( j == 1 ) exit
      end do
    end do
  else
    if ( sort == 3 ) iscore = 1
    if ( sort == 4 ) iscore = 2
    do i = 1, n-1
      j = i + 1
      do while( model(j-1)%score(iscore) < model(j)%score(iscore) )
        modeltemp = model(j-1)
        model(j-1) = model(j) 
        model(j) = modeltemp
        j = j - 1
        if ( j == 1 ) exit
      end do
    end do
  end if

end subroutine sort_by_score
