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

  character(len=4) :: refname, name

  ! If the reference atom is the CA atom, skip all except the backbone
  ! atoms remaining

  if ( refname == "CA" ) then
    if ( name /= "N" .or. &
         name /= "C" .or. &
         name /= "O" ) then
      skipatom = .true.
    end if

  ! If the reference atom is the CB atom, skip CB and side chain 

  else if ( refname == "CB" ) then
    if ( name /= "N" .or. &
         name /= "CA" .or. &
         name /= "C" .or. &
         name /= "O" ) then
      skipatom = .true.
    end if

  ! If the reference is some atom from the side chain, skip the 
  ! backbone and CB atoms

  else
    if ( name /= "N" .or. &
         name /= "CA" .or. &
         name /= "C" .or. &
         name /= "O" .or. &
         name /= "CB" ) then
      skipatom = .true.
    end if
  end if

end function skipatom

