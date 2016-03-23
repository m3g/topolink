
module topolink_operations

  implicit none

  interface operator(.eq.)
    module procedure link_eq_link, atom_eq_atom, observed_eq_observed, &
                     residue_eq_residue
  end interface 

  interface operator(.in.)
    module procedure atom_in_linktype, observed_in_linktype, &
                     atom_in_observed, atom_in_deadend, atom_in_residue, &
                     residue_in_residue, atom_in_atom
  end interface

  interface operator(.matches.)
    module procedure link_matches_observedlink, link_matches_linktype, &
                     observed_matches_linktype
  end interface

  contains

    function atom_eq_atom(atom1,atom2)
      
      use topolink_data
      implicit none
      type(pdbatom), intent(in) :: atom1, atom2
      logical :: atom_eq_atom
      
      atom_eq_atom = .false.
      if ( ( atom1%name == atom2%name ) .and. &
           ( atom1%residue%name == atom2%residue%name ) .and. &
           ( atom1%residue%chain == atom2%residue%chain ) .and. &
           ( atom1%residue%index == atom2%residue%index ) ) then
        atom_eq_atom = .true.
      end if

    end function atom_eq_atom

    function residue_eq_residue(residue1,residue2)

      use topolink_data
      implicit none
      type(pdbresidue), intent(in) :: residue1, residue2
      logical :: residue_eq_residue

      residue_eq_residue = .false.
      if ( residue1%index == residue2%index .and. &
           residue1%name == residue2%name .and. &
           residue1%chain == residue2%chain ) then
        residue_eq_residue = .true.
        return
      end if

    end function residue_eq_residue

    function link_eq_link(link1,link2)
       
      use topolink_data
      implicit none
      type(specific_link), intent(in) :: link1, link2
      logical :: link_eq_link

      link_eq_link = .false.
      if ( link1%atom1 .eq. link2%atom1 .and. &
           link1%atom2 .eq. link2%atom2 ) then 
        link_eq_link = .true.
        return
      end if
      if ( link1%atom2 .eq. link2%atom1 .and. &
           link1%atom1 .eq. link2%atom2 ) then 
        link_eq_link = .true.
        return
      end if

    end function link_eq_link

    function observed_eq_observed(observed1,observed2)

      use topolink_data
      implicit none
      type(observed_link), intent(in) :: observed1, observed2
      logical :: observed_eq_observed

      observed_eq_observed = .false.
      if ( ( observed1%residue1 .eq. observed2%residue1 ) .and. &
           ( observed1%residue2 .eq. observed2%residue2 ) ) then
        observed_eq_observed = .true.
        return
      end if
      if ( ( observed1%residue1 .eq. observed2%residue2 ) .and. &
           ( observed1%residue2 .eq. observed2%residue1 ) ) then
        observed_eq_observed = .true.
        return
      end if
      
    end function observed_eq_observed

    function atom_in_residue(atom,residue)

      use topolink_data
      implicit none
      type(pdbatom), intent(in) :: atom
      type(pdbresidue), intent(in) :: residue
      logical :: atom_in_residue

      atom_in_residue = .false.
      if ( ( atom%residue%name == residue%name .or. residue%name == 'all' ) .and. &
           ( atom%residue%chain == residue%chain .or. residue%chain == 'all' ) .and. &
           ( atom%residue%index == residue%index .or. residue%index == -1 ) ) then
        atom_in_residue = .true.
      end if

    end function atom_in_residue

    function atom_in_atom(atom1,atom2)

      use topolink_data
      implicit none
      type(pdbatom), intent(in) :: atom1, atom2
      logical :: atom_in_atom
      
      atom_in_atom = .false.
      if ( ( atom1 .in. atom2%residue ) .and. atom1%name == atom2%name ) then
        atom_in_atom = .true.
      end if

    end function atom_in_atom

    function residue_in_residue(residue1,residue2)

      use topolink_data
      implicit none
      type(pdbresidue), intent(in) :: residue1, residue2
      logical :: residue_in_residue
      
      residue_in_residue = .false.
      if ( ( residue1%name == residue2%name .or. residue2%name == 'all' ) .and. &
           ( residue1%chain == residue2%chain .or. residue2%chain == 'all' ) .and. &
           ( residue1%index == residue2%index .or. residue2%index == -1 ) ) then
        residue_in_residue = .true. 
      end if

    end function residue_in_residue

    function atom_in_linktype(atom,linktype)

      use topolink_data
      implicit none
      type(pdbatom), intent(in) :: atom
      type(link_type), intent(in) :: linktype
      logical :: atom_in_linktype

      atom_in_linktype = .false.
      if ( ( atom .in. linktype%atom1 ) .or. &
           ( atom .in. linktype%atom2 ) ) then
        atom_in_linktype = .true.
      end if

    end function atom_in_linktype

    function residue_in_linktype(residue,linktype)
      
      use topolink_data
      implicit none
      type(pdbresidue), intent(in) :: residue
      type(link_type), intent(in) :: linktype
      logical :: residue_in_linktype

      residue_in_linktype = .false.
      if ( ( residue .in. linktype%atom1%residue ) .or. &
           ( residue .in. linktype%atom2%residue ) ) then
        residue_in_linktype = .true.
      end if

    end function residue_in_linktype
  
    function observed_in_linktype(observed,linktype)

      use topolink_data
      implicit none
      type(link_type), intent(in) :: linktype
      type(observed_link), intent(in) :: observed
      logical :: observed_in_linktype

      observed_in_linktype = .false.
      if ( ( ( observed%residue1 .in. linktype%atom1%residue ) .and. &
             ( observed%residue2 .in. linktype%atom2%residue ) ) .or. &
           ( ( observed%residue1 .in. linktype%atom2%residue )   .and. &
             ( observed%residue2 .in. linktype%atom1%residue ) ) ) then
         observed_in_linktype = .true.
      end if
              
    end function observed_in_linktype

    function atom_in_observed(atom,observed)
      
      use topolink_data
      implicit none
      type(pdbatom), intent(in) :: atom
      type(observed_link), intent(in) :: observed
      logical :: atom_in_observed
      
      atom_in_observed = .false.
      if ( ( atom .in. observed%residue1 ) .or. &
           ( atom .in. observed%residue2 ) ) then
        atom_in_observed = .true.
      end if

    end function atom_in_observed

    function atom_in_deadend(atom,deadend)

      use topolink_data
      implicit none
      type(pdbatom), intent(in) :: atom
      type(observed_deadend), intent(in) :: deadend
      logical :: atom_in_deadend

      atom_in_deadend = .false.
      if ( atom .in. deadend%residue ) then 
        atom_in_deadend = .true.
      end if

    end function atom_in_deadend

    function link_matches_observedlink(link,observed)

      use topolink_data
      implicit none
      type(specific_link), intent(in) :: link
      type(observed_link), intent(in) :: observed
      logical :: link_matches_observedlink

      link_matches_observedlink = .false.
      if ( ( link%atom1 .in. observed%residue1 ) .and. &
           ( link%atom2 .in. observed%residue2 ) ) then
        link_matches_observedlink = .true.
        return
      end if
      if ( ( link%atom2 .in. observed%residue1 ) .and. &
           ( link%atom1 .in. observed%residue2 ) ) then
        link_matches_observedlink = .true.
        return
      end if

    end function link_matches_observedlink

    function link_matches_linktype(link,linktype)

      use topolink_data
      implicit none
      type(specific_link), intent(in) :: link
      type(link_type), intent(in) :: linktype
      logical :: link_matches_linktype

      link_matches_linktype = .false.
      if ( ( link%atom1 .in. linktype%atom1 ) .and. &
           ( link%atom2 .in. linktype%atom2 ) ) then
        link_matches_linktype = .true.
        return    
      end if
      if ( ( link%atom2 .in. linktype%atom1 ) .and. &
           ( link%atom1 .in. linktype%atom2 ) ) then
        link_matches_linktype = .true.
        return    
      end if

    end function link_matches_linktype

    function observed_matches_linktype(observed,linktype)

      use topolink_data
      implicit none
      type(observed_link), intent(in) :: observed
      type(link_type), intent(in) :: linktype
      logical :: observed_matches_linktype

      observed_matches_linktype = .false.
      if ( ( ( observed%residue1 .in. linktype%atom1%residue ) .and. &
             ( observed%residue2 .in. linktype%atom2%residue ) ) .or. &
             ( observed%residue2 .in. linktype%atom1%residue ) .and. &
             ( observed%residue1 .in. linktype%atom2%residue ) ) then
        observed_matches_linktype = .true.     
      end if

    end function observed_matches_linktype

end module topolink_operations











