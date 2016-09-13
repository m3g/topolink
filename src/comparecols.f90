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

  type data
    character(len=200) :: name
    double precision :: score1, score2
    logical :: foundpair
  end type data

end module comparecols_data

program comparecols

  use ioformat 
  use string_operations
  use comparecols_data
  implicit none
  integer :: imodel, i, ioerr
  integer :: nargs, nmodels, scol1, scol2, ncol1, ncol2, sort
  double precision :: score
  character(len=200) :: file1, file2, output, record, string, name
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
    write(*,*) ' -s1: The output will be sorted according to the score of the first file. ' 
    write(*,*) ' -s2: The output will be sorted according to the score of the second file. ' 
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
    write(*,*) ' ERROR: Could not read sort parameter -s0, -s1, -s2 '
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
    read(10,"(a200)",iostat=ioerr) record
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
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) (string, i = 1, scol1)
    if ( ioerr /= 0 ) cycle
    read(string,*,iostat=ioerr) score
    if ( ioerr /= 0 ) cycle
    read(record,*,iostat=ioerr) (name, i = 1, ncol1)
    if ( ioerr /= 0 ) cycle
    imodel = imodel + 1
    model(imodel)%name = name
    model(imodel)%score1 = score
  end do
  close(10)

  ! Reading data from file 2

  open(10,file=file2,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(file1))
    stop
  end if
  do imodel = 1, nmodels
    model(imodel)%foundpair = .false.
  end do
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) (string, i = 1, scol2)
    if ( ioerr /= 0 ) cycle
    read(string,*,iostat=ioerr) score
    if ( ioerr /= 0 ) cycle
    read(record,*,iostat=ioerr) (name, i = 1, ncol2)
    if ( ioerr /= 0 ) cycle
    name = basename(name)
    imodel = model_index(name,model,nmodels,error)
    if ( .not. error ) then
      model(imodel)%foundpair = .true.
      model(imodel)%score2 = score
    end if
  end do
  close(10)

  open(10,file=output)
  do imodel = 1, nmodels
    if ( model(imodel)%foundpair ) then
      write(10,"( f12.5, tr2, f12.5, tr2, a )") model(imodel)%score1, model(imodel)%score2,&
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
 
  use comparecols_data
  implicit none
  integer :: model_index
  integer :: n, imax, imin, iavg
  character(len=200) :: name
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














