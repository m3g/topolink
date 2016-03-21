!
! Module that carries the parameters needed for distance computations
! using linked cell method
!

module linkedcells

  integer :: nboxesx, nboxesy, nboxesz
  double precision :: xmin, ymin, zmin, xmax, ymax, zmax
  integer, allocatable :: ifirstbox(:,:,:), inextbox(:)

end module linkedcells
