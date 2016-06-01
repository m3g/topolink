!  
!  Written by Leandro Martínez, 2009-2011.
!  Copyright (c) 2009-2011, Leandro Martínez, Jose Mario Martinez,
!  Ernesto G. Birgin.
!  
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!  

!
! Arrays required by the flashsort package. Used only in heuristics, but
! defined here to be allocated dynamically
!

module flashsort

  implicit none
  integer, allocatable :: indflash(:) ! (ndata
  integer, allocatable :: lflash(:) ! (ndata)
  integer :: mflash ! mflash = 1 + ndata/10

end module flashsort

