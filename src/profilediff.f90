!
! Program profilediff: Computes the difference between two curves
!                      provied by discrete sets of data
! 
! L. Martinez
! Institute of Chemistry - University of Campinas
! Aug 5, 2016
!
! Run with:  ./profilediff file1.dat file2.dat
!
! Where data files are plain x y ascii tables.
!
! For greater interpolation and integral precision, increase gridsize
!

! Module that defines the maximum number of interpolation points

module size
  implicit none
  integer, parameter :: gridsize = 10000
end module size

! Module that defines the input data type

module input
  implicit none
  type input_data
    character(len=200) :: file
    double precision :: xmin, xmax
  end type input_data
end module input

!
! Main program
!

program profilediff

  use size
  use input
  implicit none
  integer :: narg, i
  double precision :: xmin, xmax, step
  double precision :: f1(gridsize), f2(gridsize), g(gridsize)
  double precision :: integral, minx
  type(input_data) :: data1, data2

  ! Read file names

  narg = iargc()
  if ( narg /= 2 ) then
    write(*,*) ' ERROR: Run with ./profilediff file1.dat file2.dat '
    stop
  end if
  call getarg(1,data1%file)
  call getarg(2,data2%file)

  ! Find xmin and xmax  

  call xminxmax(data1)
  write(*,"( 3a,f12.5 )") '# Minimum value in file ', trim(adjustl(data1%file)), ' = ', data1%xmin
  write(*,"( 3a,f12.5 )") '# Maximum value in file ', trim(adjustl(data1%file)), ' = ', data1%xmax
  call xminxmax(data2)
  write(*,"( 3a,f12.5 )") '# Minimum value in file ', trim(adjustl(data2%file)), ' = ', data2%xmin
  write(*,"( 3a,f12.5 )") '# Maximum value in file ', trim(adjustl(data2%file)), ' = ', data2%xmax

  xmin = min(data1%xmin,data2%xmin)
  xmax = max(data1%xmax,data2%xmax)
  write(*,"( a, f12.5 )") '# xmin = ', xmin
  write(*,"( a, f12.5 )") '# xmax = ', xmax
  step = ( xmax - xmin ) / dble(gridsize)
  write(*,"( a, f12.5 )") '# step = ', step
  xmin = xmin - step
  xmax = xmax + step

  call interpolatef(f1,step,xmin,data1)
  call interpolatef(f2,step,xmin,data2)

  ! Compute profilediff

  do i = 1, gridsize
    g(i) = f2(i) - f1(i)
  end do

  ! Integrate the profile difference for scores greater than 30.

  minx = 30.d0
  write(*,"( a, f8.3, a, e12.5 )") "# Integral for score greater than ", minx, ": ", &
                                   integral(g,step,xmin,minx)

  ! Write the profile difference
  
  do i = 1, gridsize 
    write(*,*) xmin+(i-1)*step, g(i)
  end do

end program profilediff

! 
! Function that integrates the functions
!

double precision function integral(f,step,xmin,minx)

  use size
  implicit none
  integer :: i
  double precision :: f(gridsize), step, xmin, minx, xgrid

  integral = 0.d0
  do i = 2, gridsize
    xgrid = xmin + (i-1)*step
    if ( xgrid < minx ) cycle
    integral = integral + step*((f(i)+f(i-1))/2.d0)
  end do

end function integral

!
! Subroutine that interpolates the input data
!

subroutine interpolatef(f,step,xmin,data)

  use size
  use input
  implicit none
  integer :: ioerr, i
  double precision :: f(gridsize), xread, y, x, slope, step, xmin
  character(len=10) :: record
  type(input_data) :: data

  open(10,file=data%file,status="old",action="read",iostat=ioerr)
  i = 1
  f(i) = 0.d0
  x = xmin
  do 
    read(10,"(a10)",iostat=ioerr) record
    if( ioerr /= 0 ) exit
    read(10,*,iostat=ioerr) xread, y
    if ( ioerr /= 0 ) cycle
    do while( x + step < xread .and. i < gridsize ) 
      slope = ( y - f(i) ) / ( xread - x ) 
      i = i + 1
      if ( x < data%xmin .or. x > data%xmax ) then
        ! Put zeros if there is no data read for this x
        f(i) = 0.d0
      else
        ! Otherwise, interpolate
        f(i) = f(i-1) + slope * step
      end if
      x = x + step
    end do
  end do
  close(10)

end subroutine interpolatef

!
! Subroutine that determines the maximum and minimum x values of the data
!

subroutine xminxmax(data)

  use input
  implicit none
  integer :: ioerr
  character(len=10) :: record
  double precision :: x
  type(input_data) :: data

  data%xmin = 1.d30
  data%xmax = -1.d30
  open(10,file=data%file,status="old",action="read",iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(data%file))
    stop
  end if
  do
    read(10,"(a10)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    read(record,*,iostat=ioerr) x
    if ( ioerr /= 0 ) cycle
    data%xmin = dmin1(x,data%xmin)
    data%xmax = dmax1(x,data%xmax)
  end do
  close(10)
  
end subroutine xminxmax


