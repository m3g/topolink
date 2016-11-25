!
! Module that carries the parameters needed for function evaluation
!

module functionpars

  integer :: natoms, atom1, atom2, nlinkatoms
  double precision :: kbond, dbond, vdwrad, kvdw, vdwrad2, dbond2, dmin
  double precision, allocatable :: coor(:,:)
  logical, allocatable :: skip(:)

end module functionpars
