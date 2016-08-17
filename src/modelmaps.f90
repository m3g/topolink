!
!
! Run with: ./modelmaps filelist2_X.txt filelist1_Y.txt filelist2_X.txt filelist2_Y.dat out1.dat out2.dat out3.dat
!
! the filelist are files containing list of LovoAlign LOG files, containing the 
! results of the alignment of each model to the first reference (X) and to the
! second reference (Y)
!
!
module datatypes

  type input_data
    character(len=200) :: xfile_list, yfile_list, plot
    character(len=200), allocatable :: xfile(:), yfile(:)
    integer :: ndata
    double precision, allocatable :: x(:), y(:)
  end type input_data
 
end module datatypes

module scoretype

  integer, parameter :: iscoretype = 1

end module scoretype

program modelmaps

  use datatypes
  use scoretype
  implicit none
  integer :: nx, ny, i, j, ix, iy
  double precision :: xmin, xmax, ymin, ymax, xgrid, ygrid
  double precision :: gridsize, zdiff_max, dist, xzmax, yzmax
  double precision, allocatable :: z1(:,:), z2(:,:), zdiff(:,:)
  character(len=200) :: diffplot
  type(input_data) :: data1, data2

  !
  ! Open files and check number of data points
  !

  call getarg(1,data1%xfile_list)
  call getndata(data1%xfile_list,data1%ndata)
  call getarg(2,data1%yfile_list)
  call getndata(data1%yfile_list,ny)
  if ( ny /= data1%ndata ) then
    write(*,*) ' ERROR: X and Y with different number of data points. '
    stop
  end if

  call getarg(3,data2%xfile_list)
  call getndata(data2%xfile_list,data2%ndata)
  call getarg(4,data2%yfile_list)
  call getndata(data2%yfile_list,ny)
  if ( ny /= data2%ndata ) then
    write(*,*) ' ERROR: X and Y with different number of data points. '
    stop
  end if

  allocate(data1%xfile(data1%ndata),data1%yfile(data1%ndata),&
           data1%x(data1%ndata),data1%y(data1%ndata),&
           data2%xfile(data2%ndata),data2%yfile(data2%ndata),&
           data2%x(data2%ndata),data2%y(data2%ndata))

  !
  ! Names of output files
  !

  call getarg(5,data1%plot)
  call getarg(6,data2%plot)
  call getarg(7,diffplot)

  !
  ! Read data
  !

  call readdata(data1)
  call readdata(data2)

  ! Remove self-alignment from data sets

  call removeself(data1%ndata,data1%x)
  call removeself(data1%ndata,data1%y)
  call removeself(data2%ndata,data2%x)
  call removeself(data2%ndata,data2%y)

  !
  ! Compute xmax and xmin for all data, to define the grid
  !

  xmin = 1.d30
  xmax = -1.d30
  ymin = 1.d30
  ymax = -1.d30
  do i = 1, data1%ndata
    xmin = dmin1(data1%x(i),xmin)
    xmax = dmax1(data1%x(i),xmax)
    ymin = dmin1(data1%y(i),ymin)
    ymax = dmax1(data1%y(i),ymax)
  end do
  do i = 1, data2%ndata
    xmin = dmin1(data2%x(i),xmin)
    xmax = dmax1(data2%x(i),xmax)
    ymin = dmin1(data2%y(i),ymin)
    ymax = dmax1(data2%y(i),ymax)
  end do

  !
  ! Number of grid points in each direction
  !
  
  gridsize = 1.d0
  nx = int((xmax-xmin)/gridsize) + 1 
  ny = int((ymax-ymin)/gridsize) + 1
  allocate(z1(nx,ny), z2(nx,ny), zdiff(nx,ny))

  ! 
  ! Populate grids
  !

  do i = 1, nx
    do j = 1, ny
      z1(i,j) = 0.d0 
      z2(i,j) = 0.d0 
      zdiff(i,j) = 0.d0 
    end do
  end do

  do i = 1, data1%ndata
    ix = int((data1%x(i)-xmin)/gridsize) + 1
    iy = int((data1%y(i)-ymin)/gridsize) + 1
    z1(ix,iy) = z1(ix,iy) + 1.d0 / data1%ndata
  end do
  
  do i = 1, data2%ndata
    ix = int((data2%x(i)-xmin)/gridsize) + 1
    iy = int((data2%y(i)-ymin)/gridsize) + 1
    z2(ix,iy) = z2(ix,iy) + 1.d0 / data2%ndata
  end do

  !
  ! Normalize differences
  !

  zdiff_max = 0.d0
  do i = 1, nx
    do j = 1, ny
      zdiff(i,j) = z2(i,j) - z1(i,j)
      zdiff_max = dmax1(zdiff_max,dabs(zdiff(i,j)))
    end do
  end do

  ! Print first data set grid
  
  open(10,file=data1%plot)
  do i = 1, nx
    do j = 1, ny
      xgrid = xmin + (i-1)*gridsize + 0.5d0*gridsize
      ygrid = ymin + (j-1)*gridsize + 0.5d0*gridsize
      if ( dabs(z1(i,j)) > 0.d0 ) then
        write(10,*) xgrid, ygrid, z1(i,j)
      end if
    end do
  end do
  close(10)

  ! Print second data set grid
  
  open(10,file=data2%plot)
  do i = 1, nx
    do j = 1, ny
      xgrid = xmin + (i-1)*gridsize + 0.5d0*gridsize
      ygrid = ymin + (j-1)*gridsize + 0.5d0*gridsize
      if ( dabs(z2(i,j)) > 0.d0 ) then
        write(10,*) xgrid, ygrid, z2(i,j)
      end if
    end do
  end do
  close(10)

  ! Printing difference grid
  
  open(10,file=diffplot)
  do i = 1, nx
    do j = 1, ny
      zdiff(i,j) = zdiff(i,j) / zdiff_max
      xgrid = xmin + (i-1)*gridsize + 0.5d0*gridsize
      ygrid = ymin + (j-1)*gridsize + 0.5d0*gridsize
      if ( dabs(zdiff(i,j)) > 0.d0 ) then
        write(10,*) xgrid, ygrid, zdiff(i,j)
      end if
    end do
  end do
  close(10)

  !
  ! Find models close to differential maximum
  !

  zdiff_max = -1.d30
  do i = 1, nx
    do j = 1, ny
      if ( zdiff(i,j) > zdiff_max ) then
        zdiff_max = zdiff(i,j)
        xgrid = xmin + (i-1)*gridsize + 0.5d0*gridsize
        ygrid = ymin + (j-1)*gridsize + 0.5d0*gridsize
        xzmax = xgrid
        yzmax = ygrid
      end if
    end do
  end do
  do i = 1, data2%ndata
    dist = (data2%x(i)-xzmax)**2 + (data2%y(i)-yzmax)**2
    dist = dsqrt(dist)
    if ( dist < 0.5d0 ) then
      write(*,*) ' Model: ', trim(adjustl(data2%xfile(i))), ' DIST: ', dist
    end if
  end do
  write(*,*) xzmax, yzmax

end program modelmaps

!
! Subroutine that reads the data points
!

subroutine readdata(data)

  use scoretype
  use datatypes
  implicit none

  integer :: i, ioerr
  character(len=200) :: record
  type(input_data) :: data

  character(len=6) :: match
  integer :: imatch, jmatch, iread, jread

  if ( iscoretype == 1 ) then
    match = "GDT_TM"
    imatch = 3
    jmatch = 8
    iread = 17
    jread = 24
  end if
  if ( iscoretype == 2 ) then
    match = "GDT_HA"
    imatch = 35
    jmatch = 40
    iread = 49
    jread = 56
  end if
  if ( iscoretype == 3 ) then
    match = "RMSD: "
    imatch = 47
    jmatch = 52
    iread = 53
    jread = 62
  end if

  open(10,file=data%xfile_list,status='old',action='read')
  i = 0 
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    i = i + 1
    data%xfile(i) = record
    open(20,file=data%xfile(i),status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not find file: ', trim(adjustl(data%xfile(i)))
      stop
    end if
    do
      read(20,"(a200)") record
      if ( record(imatch:jmatch) == match ) then
        read(record(iread:jread),*) data%x(i)
        exit
      end if
    end do
    close(20)
  end do
  close(10)
  
  open(10,file=data%yfile_list,status='old',action='read')
  i = 0 
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    i = i + 1
    data%yfile(i) = record
    open(20,file=data%yfile(i),status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not find file: ', trim(adjustl(data%yfile(i)))
      stop
    end if
    do
      read(20,"(a200)") record
      if ( record(imatch:jmatch) == match ) then
        read(record(iread:jread),*) data%y(i)
        exit
      end if
    end do
    close(20)
  end do
  close(10)

end subroutine readdata

!
! Subroutine that reads the number of data points per set
!

subroutine getndata(file,ndata)

  implicit none
  integer :: ndata, ioerr
  character(len=200) :: file, record

  open(10,file=file,status='old',action='read',iostat=ioerr) 
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(file))
  end if
  ndata = 0
  do
    read(10,*,iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    ndata = ndata + 1
  end do
  close(10)
  if ( ndata == 0 ) then
    write(*,*) ' ERROR: Could not read data from file: ', trim(adjustl(file))
    stop
  end if

end subroutine getndata

!
! Subroutine that removes self-alignment from data
!

subroutine removeself(ndata,x)

  use scoretype
  implicit none
  integer :: ndata, i, j
  double precision :: x(ndata), xtemp

  ! For scores with maximum value of 100.d0

  if ( iscoretype == 1 .or. iscoretype == 2 ) then
    j = ndata
    i = 1
    do while(i < j) 
      i = i + 1
      do while(x(i) > 99.99d0) 
        xtemp = x(i)
        x(i) = x(j)
        x(j) = xtemp
        j = j - 1
      end do
    end do
    ndata = j
  end if

  ! For RMSD

  if ( iscoretype == 3 ) then
    j = ndata
    i = 1
    do while(i < j) 
      i = i + 1
      do while(x(i) < 1.d-2) 
        xtemp = x(i)
        x(i) = x(j)
        x(j) = xtemp
        j = j - 1
      end do
    end do
    ndata = j
  end if

end subroutine removeself

  


