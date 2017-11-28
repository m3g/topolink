!
! subroutine setsigma: sets the linker vdw parameters
!
! L. Martinez
! Institute of Chemistry - University of Campinas
! http://leandro.iqm.unicamp.br
! Nov 26, 2016
!
subroutine setsigma(link,mimicchain)

  use topolink_data
  use functionpars
  implicit none
  integer :: i
  integer :: length1, length2, chainlength
  logical :: mimicchain
  type(specific_link) :: link

  ! Default radii of linker atoms is solvent radii

  do i = 1, nlinkatoms
    sigma(i) = vdwrad
  end do
  if ( .not. mimicchain ) return

  ! For CA atoms

  length1 = chainlength(link%atom1%residue%name )
  length2 = chainlength(link%atom2%residue%name )

  ! For CB atoms

  if ( link%atom1%name == "CB" ) length1 = length1 - 1
  if ( link%atom2%name == "CB" ) length2 = length2 - 1

  ! For N atoms

  if ( link%atom1%name == "N" ) length1 = length1 + 1
  if ( link%atom2%name == "N" ) length2 = length2 + 1

  ! For C atoms

  if ( link%atom1%name == "C" ) length1 = length1 + 1
  if ( link%atom2%name == "C" ) length2 = length2 + 1

  ! For O atoms

  if ( link%atom1%name == "O" ) length1 = length1 + 2
  if ( link%atom2%name == "O" ) length2 = length2 + 2

  ! For side chain atoms

  if ( link%atom1%name /= "CA" .and. &
       link%atom1%name /= "CB" .and. &
       link%atom1%name /= "N"  .and. &
       link%atom1%name /= "C"  .and. &
       link%atom1%name /= "O" ) then
    length1 = 0
  end if
  if ( link%atom2%name /= "CA" .and. &
       link%atom2%name /= "CB" .and. &
       link%atom2%name /= "N"  .and. &
       link%atom2%name /= "C"  .and. &
       link%atom2%name /= "O"  ) then
    length2 = 0
  end if

  ! Mimic atomic radii at begining of linker
write(*,*) length1, length2
  do i = 1, length1
    sigma(i) = 2.d0
  end do
  ! Mimic atomic radii at end of linker
  do i = nlinkatoms, nlinkatoms - length2 + 1, -1
    sigma(i) = 2.d0
  end do

end subroutine setsigma

integer function chainlength(name)

  implicit none
  character(len=4) :: name

  select case ( name ) 
    case ( "ALA" ) ; chainlength = 1
    case ( "ARG" ) ; chainlength = 6
    case ( "ASN" ) ; chainlength = 3
    case ( "ASP" ) ; chainlength = 4
    case ( "CYS" ) ; chainlength = 2
    case ( "GLU" ) ; chainlength = 4
    case ( "GLN" ) ; chainlength = 4
    case ( "GLY" ) ; chainlength = 1 ! Obs: needed for linker to not overlap with C and N atoms
    case ( "HIS" ) ; chainlength = 3
    case ( "HSE" ) ; chainlength = 3
    case ( "HSD" ) ; chainlength = 3
    case ( "ILE" ) ; chainlength = 3
    case ( "LEU" ) ; chainlength = 3
    case ( "LYS" ) ; chainlength = 5
    case ( "MET" ) ; chainlength = 4
    case ( "PHE" ) ; chainlength = 4
    case ( "PRO" ) ; chainlength = 2
    case ( "SER" ) ; chainlength = 2
    case ( "THR" ) ; chainlength = 2
    case ( "TRP" ) ; chainlength = 5
    case ( "TYR" ) ; chainlength = 6
    case ( "VAL" ) ; chainlength = 2
    case default ; chainlength = 0
  end select

end function chainlength
