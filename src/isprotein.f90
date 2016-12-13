!
! function isprotein: check if an atom is a protein atom
!
! L. Martinez
! Institute of Chemistry - University of Campinas
! http://leandro.iqm.unicamp.br
! Dec 13, 2016
!
logical function isprotein(atom)

  use topolink_data
  implicit none
  type(pdbatom) :: atom

  isprotein = .false.
  select case ( atom%residue%name ) 
    case ( "ALA" ) ; isprotein = .true. ; return
    case ( "ARG" ) ; isprotein = .true. ; return
    case ( "ASN" ) ; isprotein = .true. ; return
    case ( "ASP" ) ; isprotein = .true. ; return
    case ( "ASX" ) ; isprotein = .true. ; return
    case ( "CYS" ) ; isprotein = .true. ; return
    case ( "GLU" ) ; isprotein = .true. ; return
    case ( "GLN" ) ; isprotein = .true. ; return
    case ( "GLX" ) ; isprotein = .true. ; return
    case ( "GLY" ) ; isprotein = .true. ; return
    case ( "HIS" ) ; isprotein = .true. ; return
    case ( "HSE" ) ; isprotein = .true. ; return
    case ( "HSD" ) ; isprotein = .true. ; return
    case ( "ILE" ) ; isprotein = .true. ; return
    case ( "LEU" ) ; isprotein = .true. ; return
    case ( "LYS" ) ; isprotein = .true. ; return
    case ( "MET" ) ; isprotein = .true. ; return
    case ( "PHE" ) ; isprotein = .true. ; return
    case ( "PRO" ) ; isprotein = .true. ; return
    case ( "SER" ) ; isprotein = .true. ; return
    case ( "THR" ) ; isprotein = .true. ; return
    case ( "TRP" ) ; isprotein = .true. ; return
    case ( "TYR" ) ; isprotein = .true. ; return
    case ( "VAL" ) ; isprotein = .true. ; return
    case default ; return
  end select

end function isprotein
