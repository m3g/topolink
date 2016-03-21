!
! Subroutine that computes the complete function value
!

subroutine computef(n,x,f)

  use functionpars
  integer :: n
  double precision :: x(*), f, overlap, stretch
  
  f = 0.
  f = f + overlap(n,x)
  f = f + stretch(n,x)

end subroutine computef
