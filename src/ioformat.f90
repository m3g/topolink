!
! Module that contains some ahestetic format definitions
!

module ioformat
   
  integer, parameter :: max_string_length = 200
  character(len=*), parameter :: string_read = "( a200 )"

  character(len=*), parameter :: dashes = "( tr2,121('-') )"
  character(len=*), parameter :: hashes = "( tr2,121('#') )"
  character(len=*), parameter :: blank = " "

  character(len=max_string_length) :: str

  contains

  !
  ! Subroutine that writes a beatiful progress bar
  !
  
  subroutine progress(current,start,end)
  
    implicit none
    integer :: i
    integer :: current, start, end

    if ( current == start ) then
      write(*,"('  Progress: ',i10,' of ', i10$)") start, end
      return
    end if
    if ( current > start .and. current < end ) then
      write(*,"(24a,$)") (achar(8),i=1,24)
      write(*,"(i10,' of ',i10,$)") current, end
    end if
    if ( current == end ) then
      write(*,"(24a,$)") (achar(8),i=1,24)
      write(*,"(i10,' of ',i10)") current, end
    end if
   
  end subroutine progress

end module ioformat

