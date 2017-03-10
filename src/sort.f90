!
! Subroutines to sort models
!

!
! Sort by names
!

subroutine sort_by_name(n,model)

  use topolink_data
  implicit none
  integer :: i, j, n
  type(modeldata) :: model(n), modeltemp

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
! Subroutine that sorts the models by a given value of the output
! columno of the evalmodels tool
!

subroutine sort_by_value(n,model,sortcol)

  use ioformat
  use topolink_data
  implicit none
  integer :: i, j, n, sortcol
  double precision :: sortvalue
  type(modeldata) :: model(n), modeltemp

  do i = 1, n-1
    j = i + 1
    call progress(i,1,n)
    do while( sortvalue(model(j-1),sortcol) < sortvalue(model(j),sortcol) )
      modeltemp = model(j-1)
      model(j-1) = model(j) 
      model(j) = modeltemp
      j = j - 1
      if ( j == 1 ) exit
    end do
  end do
  call progress(n,1,n)

end subroutine sort_by_value

double precision function sortvalue(model,sortcol)
 
  use topolink_data
  implicit none
  integer :: sortcol
  type(modeldata) :: model

  if ( sortcol == 1 ) sortvalue = model%score
  if ( sortcol == 2 ) sortvalue = dble(model%nobscons)
  if ( sortcol == 3 ) sortvalue = dble(model%ntopcons)
  if ( sortcol == 4 ) sortvalue = dble(model%ntopnot)
  if ( sortcol == 5 ) sortvalue = dble(model%nmiss)
  if ( sortcol == 6 ) sortvalue = dble(model%nminmax)
  if ( sortcol == 7 ) sortvalue = model%sumscores
  if ( sortcol == 8 ) sortvalue = model%likeli
  if ( sortcol == 9 ) sortvalue = model%degree

end function sortvalue

! Subroutine that inverts the sorting

subroutine invert_sort(n,model)

  use topolink_data
  implicit none
  integer :: n, i, j
  type(modeldata) :: model(n), modeltemp

  i = 1
  j = n
  do while(j-i >= 1)
    modeltemp = model(i)
    model(i) = model(j)
    model(j) = modeltemp
    i = i + 1
    j = j - 1
  end do

end subroutine invert_sort



