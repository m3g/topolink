!
! subroutine setskip
!
! Defines, for a given link, which atoms must be considered
! invisible in the calculations (atoms of the same residues
! the link is attaching) 
!
! L. Martinez, Institute of Chemistry, University of Campinas
! Nov 25, 2016
!
subroutine setskip(link,atom)

  use topolink_data
  use functionpars
  implicit none
  integer :: i
  logical :: skipatom
  type(specific_link) :: link
  type(pdbatom) :: atom(natoms)

  ! Set skip to false to all atoms 
 
  do i = 1, natoms
    skip(i) = .false.
  end do

  ! Set skip to true depending on which atom is the reference
  ! atom of the link

  ! First link residue

  do i = link%atom1%residue%firstatom, link%atom1%residue%lastatom
    skip(i) = skipatom(link%atom1%name,atom(i)%name)
  end do
                
  ! Second link residue

  do i = link%atom2%residue%firstatom, link%atom2%residue%lastatom
    skip(i) = skipatom(link%atom2%name,atom(i)%name)
  end do
                
end subroutine setskip

!
! Function skipatom: sets if this atom is to be skipped
!
! Important: This must be refined in the case the structure
! has hydrogen atoms
!

logical function skipatom(refname,name)

  implicit none
  character(len=4) :: refname, name
 
  ! First, set to skip any atom

  skipatom = .true.

  ! If the reference atom is the CA atom, skip all except the backbone
  ! atoms remaining

  if ( refname == "CA" ) then
    if ( name == "N" ) skipatom = .false.
    if ( name == "C" ) skipatom = .false.
    if ( name == "O" ) skipatom = .false.

  ! If the reference atom is the CB atom, skip CB and side chain 

  else if ( refname == "CB" ) then
    if ( name == "N"  ) skipatom = .false.
    if ( name == "CA" ) skipatom = .false.
    if ( name == "C"  ) skipatom = .false.
    if ( name == "O"  ) skipatom = .false.

  ! If the refrence atom is a backbone atom, just skip everything

  else if ( refname == "N" .or. refname == "C" .or. refname == "O" ) then
    continue

  ! If the reference is some atom from the side chain, skip the 
  ! backbone and CB atoms

  else
    if ( name == "N"  ) skipatom = .false.
    if ( name == "CA" ) skipatom = .false.
    if ( name == "C"  ) skipatom = .false.
    if ( name == "O"  ) skipatom = .false.
    if ( name == "CB" ) skipatom = .false.
  end if

end function skipatom

