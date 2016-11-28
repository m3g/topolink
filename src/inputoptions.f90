!
! Module containing the input options
!
! L. Martinez
! Insitute of Chemistry - University of Campinas
! http://leandro.iqm.unicamp.br
! Nov 28, 2016
!

module inputoptions


  integer :: print
  integer :: compute
  integer :: readatoms

  double precision :: pgood
  double precision :: pbad
  double precision :: scorecut

  character(len=200) :: pdbfile
  character(len=200) :: readlog
  character(len=200) :: linkdir

  logical :: quitgood
  logical :: printlinks
  logical :: printnotfound
  logical :: observedscores
  logical :: mimicchain
  logical :: printaccessible

end module inputoptions
