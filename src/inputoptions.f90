!
! Module containing the input options
!
! L. Martinez
! Insitute of Chemistry - University of Campinas
! http://leandro.iqm.unicamp.br
! Nov 28, 2016
!

module inputoptions

  use ioformat, only : max_string_length

  integer :: print
  integer :: compute
  integer :: readatoms

  double precision :: pgood
  double precision :: pbad
  double precision :: scorecut

  character(len=max_string_length) :: pdbfile
  character(len=max_string_length) :: readlog
  character(len=max_string_length) :: linkdir

  logical :: quitgood
  logical :: printlinks
  logical :: printallfound
  logical :: printnotfound
  logical :: observedscores
  logical :: mimicchain
  logical :: printaccessible

end module inputoptions
