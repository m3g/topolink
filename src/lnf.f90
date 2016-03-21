!
! Function that computes the natural logarithm of the 
! factorial of a number
!

double precision function lnf(n)

  implicit none
  integer :: n, i

  lnf = 0.d0
  i = n
  do while( i > 0 ) 
    lnf = lnf + dlog(dble(i))
    i = i - 1
  end do

end function lnf

!
! Function that computes, safely, n*ln(p)
!

double precision function nlnp(n,p)

  implicit none
  integer :: n
  double precision :: p

  if ( n /= 0 ) then
    nlnp = n*dlog(p)
  else
    nlnp = 0.d0  
  end if

end function nlnp


