!
! Function that computes the overlap between the link and the structure
! using the linked-cell method
!

double precision function overlap(n,x)

  use functionpars
  use linkedcells
  implicit none
  integer n, ix, iy, iz, ibox, jbox, kbox, i
  double precision :: x(*)

  ! Compute distances

  dmin = vdwrad
  overlap = 0.d0
  do i = 1, nlinkatoms

    ix = (i-1)*3 + 1
    iy = ix + 1
    iz = ix + 2
    if ( x(ix) < xmin .or. x(ix) > xmax .or. &
         x(iy) < ymin .or. x(iy) > ymax .or. &
         x(iz) < zmin .or. x(iz) > zmax ) cycle
    ibox = int( ( x(ix) - xmin ) / vdwrad ) + 1
    jbox = int( ( x(iy) - ymin ) / vdwrad ) + 1
    kbox = int( ( x(iz) - zmin ) / vdwrad ) + 1

    ! Inside box

    call overlapcell(x(ix),x(iy),x(iz),ibox,jbox,kbox,overlap)

    ! Interactions of boxes that share faces

    call overlapcell(x(ix),x(iy),x(iz),ibox+1,jbox,kbox,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox,jbox+1,kbox,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox,jbox,kbox+1,overlap)

    call overlapcell(x(ix),x(iy),x(iz),ibox-1,jbox,kbox,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox,jbox-1,kbox,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox,jbox,kbox-1,overlap)

    ! Interactions of boxes that share axes

    call overlapcell(x(ix),x(iy),x(iz),ibox+1,jbox+1,kbox,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox+1,jbox,kbox+1,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox+1,jbox-1,kbox,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox+1,jbox,kbox-1,overlap)

    call overlapcell(x(ix),x(iy),x(iz),ibox,jbox+1,kbox+1,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox,jbox+1,kbox-1,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox,jbox-1,kbox+1,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox,jbox-1,kbox-1,overlap)

    call overlapcell(x(ix),x(iy),x(iz),ibox-1,jbox+1,kbox,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox-1,jbox,kbox+1,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox-1,jbox-1,kbox,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox-1,jbox,kbox-1,overlap)

    ! Interactions of boxes that share vertices

    call overlapcell(x(ix),x(iy),x(iz),ibox+1,jbox+1,kbox+1,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox+1,jbox+1,kbox-1,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox+1,jbox-1,kbox+1,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox+1,jbox-1,kbox-1,overlap)

    call overlapcell(x(ix),x(iy),x(iz),ibox-1,jbox+1,kbox+1,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox-1,jbox+1,kbox-1,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox-1,jbox-1,kbox+1,overlap)
    call overlapcell(x(ix),x(iy),x(iz),ibox-1,jbox-1,kbox-1,overlap)

  end do

end function overlap
