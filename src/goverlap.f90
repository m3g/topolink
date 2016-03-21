!
! Subroutine that computes the gradient of the overlap function using
! linked cells
!

subroutine goverlap(n,x,g)

  use functionpars
  use linkedcells
  implicit none
  integer :: n, i, ix, iy, iz, ibox, jbox, kbox
  double precision :: x(*), g(*)

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

    call goverlapcell(n,ix,iy,iz,x,ibox,jbox,kbox,g)

    ! Interactions of boxes that share faces

    call goverlapcell(n,ix,iy,iz,x,ibox+1,jbox,kbox,g)
    call goverlapcell(n,ix,iy,iz,x,ibox,jbox+1,kbox,g)
    call goverlapcell(n,ix,iy,iz,x,ibox,jbox,kbox+1,g)

    call goverlapcell(n,ix,iy,iz,x,ibox-1,jbox,kbox,g)
    call goverlapcell(n,ix,iy,iz,x,ibox,jbox-1,kbox,g)
    call goverlapcell(n,ix,iy,iz,x,ibox,jbox,kbox-1,g)

    ! Interactions of boxes that share axes

    call goverlapcell(n,ix,iy,iz,x,ibox+1,jbox+1,kbox,g)
    call goverlapcell(n,ix,iy,iz,x,ibox+1,jbox,kbox+1,g)
    call goverlapcell(n,ix,iy,iz,x,ibox+1,jbox-1,kbox,g)
    call goverlapcell(n,ix,iy,iz,x,ibox+1,jbox,kbox-1,g)

    call goverlapcell(n,ix,iy,iz,x,ibox,jbox+1,kbox+1,g)
    call goverlapcell(n,ix,iy,iz,x,ibox,jbox+1,kbox-1,g)
    call goverlapcell(n,ix,iy,iz,x,ibox,jbox-1,kbox+1,g)
    call goverlapcell(n,ix,iy,iz,x,ibox,jbox-1,kbox-1,g)

    call goverlapcell(n,ix,iy,iz,x,ibox-1,jbox+1,kbox,g)
    call goverlapcell(n,ix,iy,iz,x,ibox-1,jbox,kbox+1,g)
    call goverlapcell(n,ix,iy,iz,x,ibox-1,jbox-1,kbox,g)
    call goverlapcell(n,ix,iy,iz,x,ibox-1,jbox,kbox-1,g)

    ! Interactions of boxes that share vertices

    call goverlapcell(n,ix,iy,iz,x,ibox+1,jbox+1,kbox+1,g)
    call goverlapcell(n,ix,iy,iz,x,ibox+1,jbox+1,kbox-1,g)
    call goverlapcell(n,ix,iy,iz,x,ibox+1,jbox-1,kbox+1,g)
    call goverlapcell(n,ix,iy,iz,x,ibox+1,jbox-1,kbox-1,g)

    call goverlapcell(n,ix,iy,iz,x,ibox-1,jbox+1,kbox+1,g)
    call goverlapcell(n,ix,iy,iz,x,ibox-1,jbox+1,kbox-1,g)
    call goverlapcell(n,ix,iy,iz,x,ibox-1,jbox-1,kbox+1,g)
    call goverlapcell(n,ix,iy,iz,x,ibox-1,jbox-1,kbox-1,g)

  end do

end subroutine goverlap
