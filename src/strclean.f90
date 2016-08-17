!
! Subroutine that removes tab characters from string
!

subroutine strclean(record)

  character(len=200) :: record
  record = trim(adjustl(record))
  do i = 1, 200
    if ( record(i:i) == achar(9) ) record(i:i) = achar(32)
  end do

end subroutine strclean


