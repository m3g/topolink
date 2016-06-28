module topolink_data

  implicit none

  !
  ! Types of data
  ! 
  
  ! PDB residue data

  type pdbresidue

    integer :: index, firstatom, lastatom
    character(len=4) :: name, chain
  
  end type pdbresidue

  ! PDB atom data

  type pdbatom

     integer :: index 
     type(pdbresidue) :: residue
     character(len=4) :: name 
     double precision :: x, y, z

  end type pdbatom

  ! Specific link data

  type experiment_in_link

     ! observed: The atom pair was observed to link in this experiment
     ! type_reactive: The atom pair is reactive, according to the link types of this experiment
     ! obs_reactive: The atom pair is reactive, according to the observed reactivity in this experiment
     ! type_consistent: The topological distance of this atom pair is consistent with 
     !                  the observed reactivity of this atom pair, according to atom types
     ! obs_reactivity: The topological distance is consistent with the observed reactivity of this
     !                 atom pair, acoording to the observed atom reactivities
      
     logical :: observed, type_reactive, obs_reactive, type_consistent, obs_consistent

  end type experiment_in_link

  type specific_link

    ! atom1: First atom of this pair
    ! atom2: Second atom of this pair
    ! nbeads: Number of beads of the linker, for optimization
    ! status: Final classification of the linker
    ! n_type_expected: Number of expected links, according to atom types
    ! n_obs_expected: Number of expected links, according to observed reactivity
    ! n_type_consistent: Number of observations consistent with type-reactivity
    ! n_obs_consistent: Number of observations consistent with observed reactivity
    ! euclidean: euclidean distance
    ! topodist: topological distance, if computed 
    ! dmaxlink: Maximum length of the linker, according to experimental linkers used
    ! dmax: Maximum distance consistent with observations
    ! dmin: Minimum distance consistent with observations
    ! observed: Linker was observed to form in some experiment
    ! type_reactive: Both atoms are reactive according to some linker types
    ! obs_reactive: Both atoms are reactive according to some observed reactivity
    ! found: A topological distance was found
    ! exp: Data for the experiments, concerning this atom pair

    type(pdbatom) :: atom1, atom2
    integer :: nbeads, status
    integer :: n_type_expected, n_obs_expected, n_type_consistent, n_obs_consistent
    double precision :: euclidean, topodist, dmaxlink, dmax, dmin
    logical :: observed, type_reactive, obs_reactive, found
    type(experiment_in_link), allocatable :: exp(:)

  end type specific_link

  ! The data type for a link type

  type link_type
   
    type(pdbatom) :: atom1, atom2
    double precision :: dist

  end type link_type
 
  ! Observed link type

  type observed_link
      
    integer :: type
    double precision :: score
    type(pdbresidue) :: residue1, residue2

  end type observed_link

  ! Dead Ends

  type observed_deadend

    type(pdbresidue) :: residue   

  end type observed_deadend

  ! Experiment data

  type experiment_data

    integer :: nobs, ntypes, ndeadends, ngood, nbad
    integer :: nreactive_type, nreactive_obs
    integer :: nreach_type, nreach_obs
    integer :: noutreach_type, noutreach_obs
    integer :: nmiss_type, nmiss_obs
    character(len=200) :: name
    type(observed_deadend), allocatable :: deadend(:)
    type(observed_link), allocatable :: observed(:)
    type(link_type), allocatable :: linktype(:)
    double precision :: likelyhood, userlikelyhood, pgood, pbad, score

  end type experiment_data

  ! 
  !  Functions 
  ! 

  contains

     ! subroutine that reads atom information from a PDB line

     function read_atom(record,error)

       integer :: ioerr
       logical :: error
       type(pdbatom) :: read_atom
       character(len=200) :: record
       
       error = .false.
       if ( record(1:4) == "ATOM" .or. record(1:6) == "HETATM" ) then
         read(record(13:16),*,iostat=ioerr) read_atom%name
         if ( ioerr /= 0 ) error = .true.
         read(record(17:21),*,iostat=ioerr) read_atom%residue%name
         if ( ioerr /= 0 ) error = .true.
         read(record(22:22),*,iostat=ioerr) read_atom%residue%chain
         if ( ioerr /= 0 ) error = .true.
         read(record(23:26),*,iostat=ioerr) read_atom%residue%index
         if ( ioerr /= 0 ) error = .true.
         read(record(31:38),*,iostat=ioerr) read_atom%x
         if ( ioerr /= 0 ) error = .true.
         read(record(39:46),*,iostat=ioerr) read_atom%y
         if ( ioerr /= 0 ) error = .true.
         read(record(47:54),*,iostat=ioerr) read_atom%z
         if ( ioerr /= 0 ) error = .true.
       else
         error = .true.
       end if

     end function read_atom

     ! Prints the data of an atom

     function print_atom(atom)

        type(pdbatom) :: atom
        character(len=25) print_atom

        write(print_atom,"( 2(tr2,a4),tr2,i5,tr2,a4 )") &
                                atom%residue%name, &
                                atom%residue%chain, &
                                atom%residue%index, &
                                atom%name


     end function print_atom

     ! Reads the data of a computed link from a previous log file

     function read_link(record)

        implicit none
        integer :: i
        character(len=200) :: record
        character(len=13) :: charstat
        character(len=3) :: charobs
        character(len=9) :: chardmax
        type(specific_link) :: read_link
      
        read(record(9:12),*) read_link%atom1%residue%name
        read(record(14:14),*) read_link%atom1%residue%chain
        read(record(16:19),*) read_link%atom1%residue%index
        read(record(21:24),*) read_link%atom1%name
        
        read(record(26:29),*) read_link%atom2%residue%name
        read(record(31:31),*) read_link%atom2%residue%chain
        read(record(33:36),*) read_link%atom2%residue%index
        read(record(38:41),*) read_link%atom2%name
      
        read(record(43:50),*) read_link%euclidean

        read(record(64:66),*) charobs
        if ( charobs == "YES" ) then
          read_link%observed = .true.
        else
          read_link%observed = .false.
        end if

        read(record(79:87),*) chardmax
        do i = 1, 9
          if ( chardmax(i:i) == ">" ) chardmax(i:i) = " "
        end do
        read(chardmax,*) read_link%dmax

        read(record(89:101),"( a13 )") charstat

        if ( charstat == "    OK: FOUND" ) then
          read_link%status = 0
          read_link%found = .true.
          read(record(52:60),*) read_link%topodist
        end if
        if ( charstat == "   BAD: SHORT" ) then
          read_link%status = 1 
          read_link%found = .true.
          read(record(52:60),*) read_link%topodist
        end if
        if ( charstat == "    BAD: LONG" ) then
          read_link%status = 2
          read_link%found = .true.
          read(record(52:60),*) read_link%topodist
        end if
        if ( charstat == "    BAD: EUCL" ) then
          read_link%status = 3
          read_link%found = .false.
          read_link%topodist = -1.d0
        end if
        if ( charstat == "BAD: NOTFOUND" ) then
          read_link%status = 4
          read_link%found = .false.
          read_link%topodist = -1.d0
        end if
        if ( charstat == " BAD: MISSING" ) then
          read_link%status = 5
          read_link%found = .true.
          read(record(52:60),*) read_link%topodist
        end if
        if ( charstat == "     OK: LONG" ) then
          read_link%status = 6
          read_link%found = .true.
          read(record(52:60),*) read_link%topodist
        end if
        if ( charstat == "     OK: EUCL" ) then
          read_link%status = 7
          read_link%found = .false.
          read_link%topodist = -1.d0
        end if
        if ( charstat == " OK: NOTFOUND" ) then
          read_link%status = 8
          read_link%found = .false.
          read_link%topodist = -1.d0
        end if

     end function read_link

     ! Prints the data of a general link

     function print_link(link)

       type(specific_link) :: link
       character(len=90) :: print_link

       write(print_link,"(2(tr2,a4),tr2,i5,2(tr2,a4),tr2,i5)") &
                          link%atom1%residue%name, &
                          link%atom1%residue%chain, &
                          link%atom1%residue%index, &
                          link%atom2%residue%name, &
                          link%atom2%residue%chain, &
                          link%atom2%residue%index

     end function print_link

     ! Prints the data of an observed link

     function print_obs(observed)

       type(observed_link) :: observed
       character(len=90) :: print_obs

       write(print_obs,"(2(tr2,a4),tr2,i5,2(tr2,a4),tr2,i5)") &
                          observed%residue1%name, &
                          observed%residue1%chain, &
                          observed%residue1%index, &
                          observed%residue2%name, &
                          observed%residue2%chain, &
                          observed%residue2%index

     end function print_obs

     ! Prints the data of an observed dead end

     function print_deadend(deadend)

       type(observed_deadend) :: deadend
       character(len=90) :: print_deadend

       write(print_deadend,"( 2(tr2,a4),tr2,i5 )") &
                              deadend%residue%name, &
                              deadend%residue%chain, &
                              deadend%residue%index

     end function print_deadend

     ! Prints the data of a linktype

     function print_linktype(link)

       type(link_type) :: link
       character(len=90) :: print_linktype
        
       write(print_linktype,"( 6(tr2,a4),tr2,f8.3)") &
                               link%atom1%residue%name,&
                               link%atom1%residue%chain,&
                               link%atom1%name,&
                               link%atom2%residue%name,&
                               link%atom2%residue%chain,&
                               link%atom2%name,&
                               link%dist

     end function print_linktype

     ! Prints a PDB atom line with coordinates
     
     function print_pdbatom(atom)

       type(pdbatom) :: atom
       character(len=200) :: pdbformat, print_pdbatom

       pdbformat = "('ATOM',t7,i5,t12,a4,t17,a4,t22,a1,t23,i4,t31,f8.3,t39,f8.3,t47,f8.3)"
       write(print_pdbatom,pdbformat) &
                           atom%index, trim(adjustl(atom%name)), &
                           trim(adjustl(atom%residue%name)), atom%residue%chain, &
                           atom%residue%index, atom%x, atom%y, atom%z

     end function print_pdbatom

end module topolink_data

