!
! Subroutine that creates the initial guess for the path
!

subroutine initguess(n,x,iguess)

  use functionpars
  implicit none
  integer :: n, iguess, ix, iy, iz, i, j, ntrial, itrial
  double precision :: x(n), theta, phi, random, vec(3), overlap
  double precision :: current, best, xbest(n)
  double precision, parameter :: pi = 4.d0*datan(1.d0)

  double precision :: t, t_right, t_left, tstep
  double precision :: xleft(3), xright(3), xtemp(3)
  double precision :: xn(3), xnorm
  double precision :: xcenter(3), radius
  double precision :: xp1(3), p1norm
  double precision :: xp2(3), p2norm
  double precision :: t_circle, vnorm
  double precision :: dleft, dright, darch, dtot, dstep
  double precision :: vleft(3), vright(3)
  double precision :: random1, random2
  integer :: nsteps
  integer :: narch, nright, nleft
  double precision :: step

  ! Initial guesses as straight lines pointing outwards
  ! from atom1, at random directions 

  ntrial = 100

  iguess = 3
  if ( iguess == 1 ) then
    best = 1.d30
    do j = 1, ntrial
      call random_number(random)
      theta = pi*random
      call random_number(random)
      phi = 2.d0*pi*random
      vec(1) = dsin(theta)*dcos(phi) 
      vec(2) = dsin(theta)*dsin(phi)
      vec(3) = dcos(theta)
      do i = 1, nlinkatoms
        ix = (i-1)*3 + 1
        iy = ix + 1
        iz = ix + 2
        x(ix) = coor(atom1,1) + vec(1)*i*dbond
        x(iy) = coor(atom1,2) + vec(2)*i*dbond
        x(iz) = coor(atom1,3) + vec(3)*i*dbond
      end do
      current = overlap(n,x)
      if ( current < best ) then
        do i = 1, n
          xbest(i) = x(i)
        end do
      end if
      if ( overlap(n,x) < 1.d0 ) exit
    end do
    do i = 1, n
      x(i) = xbest(i) 
    end do
  end if

  ! Similar to iguess=1, but the segment is divided in two, pointing
  ! each from final and end atoms

  if ( iguess == 2 .or. iguess == 3 ) then
    best = 1.d30
    do j = 1, ntrial
      call random_number(random)
      theta = pi*random
      call random_number(random)
      phi = 2.d0*pi*random
      vec(1) = dsin(theta)*dcos(phi) 
      vec(2) = dsin(theta)*dsin(phi)
      vec(3) = dcos(theta)
      do i = 1, nlinkatoms/2
        ix = (i-1)*3 + 1
        iy = ix + 1
        iz = ix + 2
        x(ix) = coor(atom1,1) + vec(1)*i*dbond
        x(iy) = coor(atom1,2) + vec(2)*i*dbond
        x(iz) = coor(atom1,3) + vec(3)*i*dbond
      end do
      call random_number(random)
      theta = pi*random
      call random_number(random)
      phi = 2.d0*pi*random
      vec(1) = dsin(theta)*dcos(phi) 
      vec(2) = dsin(theta)*dsin(phi)
      vec(3) = dcos(theta)
      do i = nlinkatoms, nlinkatoms/2+1, -1
        ix = (i-1)*3 + 1
        iy = ix + 1
        iz = ix + 2
        x(ix) = coor(atom2,1) + vec(1)*(nlinkatoms-i+1)*dbond
        x(iy) = coor(atom2,2) + vec(2)*(nlinkatoms-i+1)*dbond
        x(iz) = coor(atom2,3) + vec(3)*(nlinkatoms-i+1)*dbond
      end do
      current = overlap(n,x)
      if ( current < best ) then
        best = current
        do i = 1, n
          xbest(i) = x(i)
        end do
      end if
      if ( current < 1.d0 ) exit
    end do
    do i = 1, n
      x(i) = xbest(i)
    end do
  end if

  ! voltar

  if ( iguess == 3 ) then

    ix = (nlinkatoms/2-1)*3+1
    xleft(1) = x(ix)
    xleft(2) = x(ix+1)
    xleft(3) = x(ix+2)
    ix = (nlinkatoms/2+1-1)*3+1
    xright(1) = x(ix)
    xright(2) = x(ix+1)
    xright(3) = x(ix+2)

    ! Compute the length of each segment, to know how many points will
    ! be in the arch

    vleft(1) = xleft(1) - coor(atom1,1)
    vleft(2) = xleft(2) - coor(atom1,2)
    vleft(3) = xleft(3) - coor(atom1,3)
    dleft = vnorm(vleft)
    vleft(1) = vleft(1) / dleft
    vleft(2) = vleft(2) / dleft
    vleft(3) = vleft(3) / dleft
    dleft = dmin1(dleft,15.d0)
    xleft(1) = coor(atom1,1) + dleft*vleft(1)
    xleft(2) = coor(atom1,2) + dleft*vleft(2)
    xleft(3) = coor(atom1,3) + dleft*vleft(3)

    vright(1) = coor(atom2,1) - xright(1)
    vright(2) = coor(atom2,2) - xright(2)
    vright(3) = coor(atom2,3) - xright(3)
    dright = vnorm(vright)
    vright(1) = vright(1) / dright
    vright(2) = vright(2) / dright
    vright(3) = vright(3) / dright
    dright = dmin1(dright,15.d0)
    xright(1) = coor(atom2,1) - dright*vright(1)
    xright(2) = coor(atom2,2) - dright*vright(2)
    xright(3) = coor(atom2,3) - dright*vright(3)

    ! xp1 is one of the vectors of the plane containing the circle

    xp1(1) = (xleft(1) - xright(1))
    xp1(2) = (xleft(2) - xright(2))
    xp1(3) = (xleft(3) - xright(3))
    p1norm = vnorm(xp1)
    xp1(1) = xp1(1) / p1norm
    xp1(2) = xp1(2) / p1norm
    xp1(3) = xp1(3) / p1norm

    ! xcenter is the center of the circle

    xcenter(1) = (xleft(1) + xright(1))/2.d0
    xcenter(2) = (xleft(2) + xright(2))/2.d0
    xcenter(3) = (xleft(3) + xright(3))/2.d0

    ! the radius is set by the distance between the endpoints
   
    radius = dsqrt((xleft(1)-xcenter(1))**2 + &
                   (xleft(2)-xcenter(2))**2 + &
                   (xleft(3)-xcenter(3))**2)

    best = 1.d30 
    do itrial = 1, ntrial

      ! Now create what will be the normal vector of the circle

      call random_number(random)
      random1 = -1.d0 + 2.d0*random
      call random_number(random)
      random2 = -1.d0 + 2.d0*random
      if ( xp1(1) == 0. .and. xp1(2) == 0. ) then ! cz = 0
        xn(1) = random1
        xn(2) = random2
        xn(3) = 0.
      else if ( xp1(1) == 0. .and. xp1(3) == 0. ) then ! by = 0
        xn(1) = random1
        xn(2) = 0.
        xn(3) = random2
      else if ( xp1(2) == 0. .and. xp1(3) == 0. ) then ! ax = 0
        xn(1) = 0.
        xn(2) = random1
        xn(3) = random2
      else if ( xp1(1) == 0. ) then  ! by + cz = 0 
        xn(1) = random1
        xn(2) = random2
        xn(3) = -xp1(2)*xn(2)/xp1(3)
      else if ( xp1(2) == 0. ) then ! ax + cz = 0
        xn(1) = random1
        xn(2) = random2
        xn(3) = -xp1(1)*xn(1)/xp1(3)
      else if ( xp1(3) == 0. ) then ! ax + by = 0
        xn(1) = random1
        xn(2) = -xp1(1)*xn(1)/xp1(2)
        xn(3) = random2
      else ! ax + by + cz = 0
        xn(1) = random1
        xn(1) = random2
        xn(3) = ( xp1(1)*xn(1) + xp1(2)*xn(2) / xp1(3) )
      end if
      xnorm = vnorm(xn)
      xn(1) = xn(1) / xnorm
      xn(2) = xn(2) / xnorm
      xn(3) = xn(3) / xnorm

      ! The second vector in the circle plane is perpendicular to
      ! both xn and xp1

      xp2(1) = ( xp1(2)*xn(3) - xp1(3)*xn(2) ) 
      xp2(2) = -1.d0*( xp1(1)*xn(3) - xp1(3)*xn(1) )
      xp2(3) = ( xp1(1)*xn(2) - xp1(2)*xn(1) )
      p2norm = vnorm(xp2)
      xp2(1) = xp2(1) / p2norm
      xp2(2) = xp2(2) / p2norm
      xp2(3) = xp2(3) / p2norm

      ! Determine which is the parametric value for xright and xleft

      t_right = t_circle(xright,xcenter,xp1,xp2)
      t_left = t_circle(xleft,xcenter,xp1,xp2)

      darch = radius * abs( t_right - t_left )
      dtot = dleft + dright + darch
      dstep = dtot / (nlinkatoms+1)

      narch = darch/dstep+1
      tstep = (t_right - t_left) / narch

      ! Build the initial guess for the linker

      i = 0
      dtot = 0.d0
      nleft = 0
      narch = 0
      do while( i < nlinkatoms ) 
        i = i + 1
        dtot = dtot + dstep
        ix = (i-1)*3 + 1
        iy = ix + 1
        iz = ix + 2
        if ( dtot <= dleft ) then
          x(ix) = coor(atom1,1) + i*dstep*vleft(1)
          x(iy) = coor(atom1,2) + i*dstep*vleft(2)
          x(iz) = coor(atom1,3) + i*dstep*vleft(3)
          nleft = nleft + 1
        else if ( dtot > dleft .and. dtot <= dleft+darch ) then
          t = t_left + (i-nleft)*tstep
          x(ix) = xcenter(1) + radius*dcos(t)*xp1(1) + radius*dsin(t)*xp2(1)
          x(iy) = xcenter(2) + radius*dcos(t)*xp1(2) + radius*dsin(t)*xp2(2)
          x(iz) = xcenter(3) + radius*dcos(t)*xp1(3) + radius*dsin(t)*xp2(3)
          narch = narch + 1
        else
          x(ix) = xright(1) + (i-nleft-narch)*dstep*vright(1)
          x(iy) = xright(2) + (i-nleft-narch)*dstep*vright(2)
          x(iz) = xright(3) + (i-nleft-narch)*dstep*vright(3)
        end if
      end do

      current = overlap(n,x)
      if ( current < best ) then
        best = current
        do i = 1, n
          xbest(i) = x(i)
        end do
      end if
      if ( current < 1.d0 ) exit

    end do
    do i = 1, n
      x(i) = xbest(i)
    end do

  end if

end subroutine initguess

!
! This function determines which is the parameter, t, of a vector
! belonging to a unitary circle parametrized as r(t) = c + a*cos(t) + b*sin(t)
! r is the vector, and "a" and "b" are ortonormal vectors in the circle plane,
! and c is the center of the cirle in 3D
!

function t_circle(r,c,a,b)

  implicit none
  double precision :: t_circle
  double precision :: r(3), c(3), a(3), b(3)
  double precision :: v(3)
  double precision :: vn
  double precision :: vtrial(3), sint, cost
  double precision :: error
  double precision :: vnorm, cramer_x
  double precision, parameter :: pi = 4.d0*datan(1.d0)

  v(1) = r(1) - c(1)
  v(2) = r(2) - c(2)
  v(3) = r(3) - c(3)
  vn = vnorm(v)
  v(1) = v(1)/vn
  v(2) = v(2)/vn
  v(3) = v(3)/vn

  if ( a(1)*b(2) - b(1)*a(2) /= 0. ) then ! eqs. 1 and 2
    cost = cramer_x(a(1),a(2),b(1),b(2),v(1),v(2))
  else if ( a(1)*b(3) - b(1)*a(3) /= 0. ) then ! eqs. 1 and 3
    cost = cramer_x(a(1),a(3),b(1),b(3),v(1),v(3))
  else if ( a(2)*b(3) - b(2)*a(3) /= 0. ) then ! eqs 2 and 3
    cost = cramer_x(a(2),a(3),b(2),b(3),v(2),v(3))
  else 
    write(*,*) ' ops, there is a problem. '
    stop
  end if
  if ( cost > 1.d0 ) cost = 1.d0
  if ( cost < -1.d0 ) cost = -1.d0
  t_circle = acos(cost) 

  ! just for checking 
  ! vtrial(1) = c(1) + vn*( cost*a(1) + dsqrt(1.d0-cost**2)*b(1) )
  ! vtrial(2) = c(2) + vn*( cost*a(2) + dsqrt(1.d0-cost**2)*b(2) )
  ! vtrial(3) = c(3) + vn*( cost*a(3) + dsqrt(1.d0-cost**2)*b(3) )
  ! error = (r(1) - vtrial(1))**2 + (r(2)-vtrial(2))**2 + (r(3)-vtrial(3))**2

end function t_circle

! Function that returns the norm of a 3D vector

function vnorm(v)

  implicit none
  double precision :: v(3), vnorm

  vnorm = dsqrt(v(1)**2 + v(2)**2 + v(3)**2)

end function vnorm

! Function that returns the value of x using Cramer's formula
! for a 2x2 linear system: [(a1,a2),(b1,b2)](x,y)=(v1,v2)

function cramer_x(a1,a2,b1,b2,c1,c2)

  implicit none
  double precision :: cramer_x
  double precision :: a1, b1, a2, b2, c1, c2

  cramer_x = ( c1*b2 - b1*c2 ) / ( a1*b2 - b1*a2 )

end function cramer_x









