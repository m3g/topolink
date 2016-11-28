!
! Subroutine that computes the complete gradient
!

subroutine computeg(n,x,g)

  use functionpars
  integer :: n, i
  double precision :: x(*), g(*)

  do i = 1, n
    g(i) = 0.d0
  end do
  call goverlap(n,x,g)
  call gstretch(n,x,g)

end subroutine computeg
