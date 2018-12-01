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
  double precision :: search_limit
  integer :: search_limit_type

  character(len=max_string_length) :: pdbfile
  character(len=max_string_length) :: readlog
  character(len=max_string_length) :: linkdir

  logical :: quitgood
  logical :: printlinks
  logical :: printallfound
  logical :: printnotfound
  logical :: printPDBnotfound
  logical :: observedscores
  logical :: mimicchain
  logical :: printaccessible

  logical :: screen_log
  logical :: output_log
  character(len=max_string_length) :: output_log_file
  integer :: output_log_unit 

end module inputoptions
