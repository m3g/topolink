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

