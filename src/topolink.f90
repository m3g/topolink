!
! Program TOPOLINK
!
! L. Martinez, 
! Institute of Chemistry, University of Campinas
! http://leandro.iqm.unicamp.br
! Out 5, 2015
!
! Reference:
!
! L. Martinez, A. Ferrari, F. C. Gozzo,
! TopoLink: A package to compute the likelyhood of structural models
! based on surface accessible topological distances
! 2015
!

program topolink

  use ioformat
  use functionpars
  use linkedcells
  use topolink_data
  use topolink_operations
  use string_operations
  use inputoptions
  implicit none
  integer :: nargs, iargc, ioerr, i, j, k, ii, n, seed, ix, iy, iz, ntrial, itrial, &
             iguess, optpars(10), best_repeat, nbest, nobs, ngooddist, nbaddist, nmisslinks, &
             i1, i2, ndeadends, nexp, iexp, ntypes, npairs, &
             ngood, natreactive, nmax, nloglines, linkstatus, first(2), last(2), nchains
  double precision:: f, stretch, overlap, dpath, dpath_best, computedpath, overviol, &
                     kpath, likelyhood, userlikelyhood, lnf, nlnp, totscore, readscore,& 
                     kvdwini
  character(len=4) :: char1, char2 
  character(len=200) :: record, linkfile, inputfile, endread
  character(len=200), allocatable :: logline(:)
  character(len=20) :: floatout, intout, intout2
  logical :: error, r1, r2, inexp, warning, interchain, isprotein

  external :: computef, computeg

  logical, allocatable :: type_reactive(:,:), obs_reactive(:,:), observed(:,:)
  integer, allocatable :: reactiveatom(:), iobserved(:,:)
  double precision, allocatable :: dlink(:,:), x(:), g(:), xbest(:)

  type(experiment_data), allocatable :: experiment(:)
  type(pdbatom) :: readatom, writeatom, testatom
  type(pdbatom), allocatable :: atom(:)
  type(specific_link), allocatable :: link(:)
  type(specific_link) :: linktest

  ! Print title

  call title()

  ! Check if the number of arguments is correct

  nargs = iargc()
  if ( nargs /= 1 .and. nargs /= 2 .and. nargs /= 3 ) then
    write(*,*) ' Run with: '
    write(*,*) 
    write(*,*) ' topolink inputfile.inp [pdbfile] '
    write(*,*) 
    write(*,hashes)
    write(*,*) 
    stop
  end if

  ! Default input parametes

  pdbfile = 'none'
  readlog = 'none' 
  print = 0
  quitgood = .false.
  linkdir = 'none'
  printlinks = .false.
  printnotfound = .false.
  compute = 2
  pgood = -1.d0
  pbad = -1.d0
  observedscores = .false.
  scorecut = -1.d0
  mimicchain = .true.
  printaccessible = .false.
  warning = .false.
  interchain = .false.

  ! Format of output
  floatout="(tr1,a,f12.7)"
  intout="(tr1,a,i8)"
  intout2="( tr1,a,i5,tr2,i5 )"

  ! Default function parameters

  dbond = 1.5
  kvdw = 2.
  kvdwini = kvdw
  kbond = 10.
  vdwrad = 3.
  ntrial = 100
  nbest = 5
  kpath = 10.d0
  endread = "##################"
  readatoms = 1
 
  ! Parameters for optimization method

  dbond2 = dbond**2
  optpars(1) = 500 ! Maximum number of functional evaluations
  optpars(2) = 100 ! Maximum number of CG iterations
  seed = 0
  iguess = 1

  ! Read command line parameters

  call getarg(1,inputfile)
  if ( nargs > 1 ) call getarg(2,pdbfile)
  if ( nargs > 2 ) call getarg(3,readlog)
  
  open(10,file=inputfile,status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open input file: ', trim(adjustl(inputfile))
    stop
  end if

  ! Read input file

  nexp = 0
  error = .false.
  input : do 
    read(10,"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    call strclean(record) 
    if ( length(record) < 1 .or. record(1:1) == "#" ) cycle

    select case ( keyword(record) )
      case ("pdbfile")
        if ( pdbfile == 'none' ) then
          pdbfile = keyvalue(record,1)
        end if
      case ("readlog")
        if ( readlog == 'none' ) then
          readlog = keyvalue(record,1)
        end if
      case ("endread")
        endread = keyvalue(record,1)
      case ("readatoms") 
        if ( keyvalue(record,1) == 'all' ) readatoms = 0
        if ( keyvalue(record,1) == 'heavy' ) readatoms = 1
        if ( keyvalue(record,1) == 'backbone' ) readatoms = 2
        if ( keyvalue(record,1) == 'backplusCB' ) readatoms = 3
      case ("linkdir")
        linkdir = keyvalue(record,1)
      case ("observed")
        cycle
      case ("deadend")
        cycle
      case ("linktype") 
        cycle
      case ("experiment")
        nexp = nexp + 1
      case ("end")
        cycle
      case ("print")
        record = keyvalue(record,1)
        read(record,*) print
      case ("mimicchain")
        if ( keyvalue(record,1) == 'yes' ) mimicchain = .true.
        if ( keyvalue(record,1) == 'no' ) mimicchain = .false.
      case ("interchain")
        interchain = .true.
      case ("quitgood")
        if ( keyvalue(record,1) == 'yes' ) quitgood = .true.
        if ( keyvalue(record,1) == 'no' ) quitgood = .false.
      case ("printlinks")
        if ( keyvalue(record,1) == 'yes' ) printlinks = .true.
        if ( keyvalue(record,1) == 'no' ) printlinks = .false.
      case ("printnotfound")
        if ( keyvalue(record,1) == 'yes' ) printnotfound = .true.
        if ( keyvalue(record,1) == 'no' ) printnotfound = .false.
      case ("printaccessible")
        if ( keyvalue(record,1) == 'yes' ) printaccessible = .true.
        if ( keyvalue(record,1) == 'no' ) printaccessible = .false.
      case ("compute")
        if ( keyvalue(record,1) == 'observed' ) compute = 1
        if ( keyvalue(record,1) == 'reactive' ) compute = 2
        if ( keyvalue(record,1) == 'all' ) compute = 3
      case ("pgood") 
        record = keyvalue(record,1)
        read(record,*) pgood
      case ("pbad") 
        record = keyvalue(record,1)
        read(record,*) pbad
      case ("dbond") 
        record = keyvalue(record,1)
        read(record,*) dbond
      case ("kbond") 
        record = keyvalue(record,1)
        read(record,*) kbond
      case ("scorecut") 
        record = keyvalue(record,1)
        read(record,*) scorecut
      case ("vdwrad") 
        record = keyvalue(record,1)
        read(record,*) vdwrad
      case ("kpath") 
        record = keyvalue(record,1)
        read(record,*) kpath
      case ("kvdw") 
        record = keyvalue(record,1)
        read(record,*) kvdw
        kvdwini = kvdw
      case ("ntrial") 
        record = keyvalue(record,1)
        read(record,*) ntrial
      case ("nbest")
        record = keyvalue(record,1)
        read(record,*) nbest
      case ("maxfunc")
        record = keyvalue(record,1)
        read(record,*) optpars(1)
      case ("maxcg")
        record = keyvalue(record,1)
        read(record,*) optpars(2)
      case ("iguess")
        record = keyvalue(record,1)
        read(record,*) iguess
      case ("seed")
        if ( keyvalue(record,1) == 'random' ) then
          seed = 0
        else
          record = keyvalue(record,1)
          read(record,*) seed
        end if
      case ("exit") 
        exit input
      case default
        error = .true.
        write(*,*) ' ERROR: Input parameter not recognized: ', keyword(record)
    end select
    
  end do input
  if ( error ) then
    close(10)
    stop
  end if
  rewind(10)

  ! Random number initialization
  
  if ( seed == 0 ) call seed_from_time(seed)
  call init_random_number(seed)
 
  ! Print input parameters

  write(*,"(a,a)") '  PDB input file: ', trim(adjustl(pdbfile))
  if ( endread /= "###" ) then
    write(*,*) ' Will stop PDB reading when ', trim(endread), ' is found. '
  end if
  if ( readatoms == 0 ) then
    write(*,*) ' All atoms will be considered, including hydrogens. '
  else if ( readatoms == 1 ) then
    write(*,*) ' Only heavy atoms will be considered. '
  else if ( readatoms == 2 ) then
    write(*,*) ' Only backbone atoms will be considered. '
  else if ( readatoms == 3 ) then
    write(*,*) ' Only backbone plus CB atoms will be considered. '
  else
    write(*,*) " ERROR: Wrong keyword value for 'readatoms' "
    stop
  end if
  write(*,*)
  if ( .not. printlinks ) then
    write(*,*) ' Links will not be written as PDB files. '
  else 
    if ( linkdir == 'none' ) then
      write(*,*) ' ERROR: Need linkdir if printlinks is yes. '
      stop
    else 
      i = length(linkdir)
      if ( linkdir(i:i) /= "/" ) then
        linkdir(i+1:i+1) = "/"
      end if
      write(*,*) ' Directory to output links as PDB files: ', trim(adjustl(linkdir))
    end if
  end if

  write(*,*) 
  write(*,intout) ' Printing option: ', print
  write(*,*) ' Leave when first valid link is found: ', quitgood
  write(*,floatout) ' Link bond distance: ', dbond
  write(*,floatout) ' Link bond force-constant: ', kbond
  write(*,floatout) ' VdW radius for volume exclusion: ', vdwrad
  write(*,floatout) ' Path energy constant: ', kpath
  write(*,floatout) ' Cutoff for observed scores: ', scorecut
  write(*,intout) ' Number of trials for link search: ', ntrial
  write(*,intout) ' Number of repeated best links to find until quit: ', nbest
  write(*,intout) ' Maximum number of function evaluations of CGNewton: ', optpars(1)
  write(*,intout) ' Maximum number of CG iterations in CGNewton: ', optpars(2)
  write(*,intout) ' Seed for random number generator: ', seed

  if ( readlog /= "none" ) then
    write(*,*) ' Will read link data from previous run: ', trim(adjustl(readlog))
    open(20,file=readlog,action='read',iostat=ioerr)
    nloglines = 0
    do      
      read(20,"( a200 )",iostat=ioerr) record
      if ( ioerr /= 0 ) exit
      call strclean(record)
      if ( record(3:7) == "LINK:" ) then
        nloglines = nloglines + 1
      end if
    end do
    rewind(20)
    allocate(logline(nloglines))
    i1 = 0
    do
      read(20,"( a200 )",iostat=ioerr) record
      if ( ioerr /= 0 ) exit
      call strclean(record)
      if ( record(3:7) == "LINK:" ) then
        i1 = i1 + 1
        logline(i1) = record
      end if
    end do
    write(*,*) ' Found ', nloglines, ' links in log file. '
    close(20)
  end if

  write(*,*)
  write(*,intout) ' Number of experiments: ', nexp
  allocate( experiment(nexp) )

  !
  ! Reading and annotating number of observed links, types of links and deadends
  ! 

  inexp = .false.
  iexp = 0
  do 
    read(10,"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    call strclean(record)
    if ( length(record) < 1 .or. record(1:1) == "#" ) cycle

    if ( keyword(record) == 'experiment' ) then
      inexp = .true.
      iexp = iexp + 1
      nobs = 0
      ntypes = 0
      ndeadends = 0
      if ( length(record) > 11 ) then
        experiment(iexp)%name = record(11:length(record))
      else
        write(record,*) iexp
        experiment(iexp)%name = record(11:length(record))
      end if
      cycle
    end if

    if ( keyword(record) == 'end' ) then
      if ( keyvalue(record,1) == 'experiment' ) then
        experiment(iexp)%nobs = nobs
        experiment(iexp)%ntypes = ntypes
        experiment(iexp)%ndeadends = ndeadends
        allocate( experiment(iexp)%observed(nobs) )
        allocate( experiment(iexp)%linktype(ntypes) )
        allocate( experiment(iexp)%deadend(ndeadends) )
        inexp = .false.
        cycle
      end if
    end if

    if ( keyword(record) == 'observed' ) then
      ! Reads only observed links with scores greater than scorecut
      readscore = 0.d0
      if ( countwords(record) == 8 ) then
        observedscores = .true.
        char1 = trim(keyvalue(record,7))
        read(char1,*) readscore
      end if
      if ( readscore < scorecut ) cycle
      nobs = nobs + 1
    end if
    if ( keyword(record) == 'linktype'  ) ntypes = ntypes + 1
    if ( keyword(record) == 'deadend' ) ndeadends = ndeadends + 1

    if ( .not. inexp ) then
      if ( keyword(record) == 'observed' .or. &
           keyword(record) == 'linktype' .or. &
           keyword(record) == 'deadend' ) then
        write(*,*) ' ERROR: On input: observed, linktype or deadend keyword outside experiment section. '
        stop
      end if
    end if

    if ( keyword(record) == 'exit' ) exit

  end do
  rewind(10)

  ! 
  ! Reading the actual data for observed links, link types and dead ends
  ! 

  iexp = 0
  do 
    read(10,"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    call strclean(record)
    if ( length(record) < 1 .or. record(1:1) == "#" ) cycle

    if ( keyword(record) == 'experiment' ) then
      iexp = iexp + 1
      nobs = 0
      ntypes = 0
      ndeadends = 0
      if ( length(record) > 11 ) then
        experiment(iexp)%name = trim(adjustl(record(11:length(record))))
      else
        write(record,*) iexp
        experiment(iexp)%name = trim(adjustl(record(11:length(record))))
      end if
      cycle
    end if

    if ( keyword(record) == 'observed' ) then 
      readscore = 0.d0
      if ( countwords(record) == 8 ) then
        observedscores = .true.
        char1 = trim(keyvalue(record,7))
        read(char1,*) readscore
      end if
      if ( readscore < scorecut ) cycle
      nobs = nobs + 1
      experiment(iexp)%observed(nobs)%residue1%name = trim(keyvalue(record,1))
      experiment(iexp)%observed(nobs)%residue1%chain = trim(keyvalue(record,2))
      char1 = trim(keyvalue(record,3))
      read(char1,*) experiment(iexp)%observed(nobs)%residue1%index
      experiment(iexp)%observed(nobs)%residue2%name = trim(keyvalue(record,4))
      experiment(iexp)%observed(nobs)%residue2%chain = trim(keyvalue(record,5))
      char1 = trim(keyvalue(record,6))
      read(char1,*) experiment(iexp)%observed(nobs)%residue2%index
      experiment(iexp)%observed(nobs)%score = readscore
    end if                      

    if ( keyword(record) == 'linktype'  ) then
      ntypes = ntypes + 1
      experiment(iexp)%linktype(ntypes)%atom1%residue%name = trim(keyvalue(record,1))
      experiment(iexp)%linktype(ntypes)%atom1%residue%chain = trim(keyvalue(record,2))
      char1 = trim(keyvalue(record,3))
      if ( char1 == 'all' ) then
        experiment(iexp)%linktype(ntypes)%atom1%residue%index = -1
      else
        read(char1,*) experiment(iexp)%linktype(ntypes)%atom1%residue%index
      end if
      experiment(iexp)%linktype(ntypes)%atom1%name = trim(keyvalue(record,4))
      experiment(iexp)%linktype(ntypes)%atom2%residue%name = trim(keyvalue(record,5))
      experiment(iexp)%linktype(ntypes)%atom2%residue%chain = trim(keyvalue(record,6))
      char1 = trim(keyvalue(record,7))
      if ( char1 == 'all' ) then
        experiment(iexp)%linktype(ntypes)%atom2%residue%index = -1
      else
        read(char1,*) experiment(iexp)%linktype(ntypes)%atom2%residue%index
      end if
      experiment(iexp)%linktype(ntypes)%atom2%name = trim(keyvalue(record,8))
      char1 = trim(keyvalue(record,9))
      read(char1,*) experiment(iexp)%linktype(ntypes)%dist
    end if

    if ( keyword(record) == 'deadend' ) then
      ndeadends = ndeadends + 1
      experiment(iexp)%deadend(ndeadends)%residue%name = trim(keyvalue(record,1))
      experiment(iexp)%deadend(ndeadends)%residue%chain = trim(keyvalue(record,2))
      char1 = trim(keyvalue(record,3))
      read(char1,*) experiment(iexp)%deadend(ndeadends)%residue%index
    end if

    if ( keyword(record) == 'end' ) then
      if( keyvalue(record,1) == 'experiment' ) cycle
    end if
    if ( keyword(record) == 'exit' ) exit

  end do
  close(10)

  !
  ! Checking if the input of observed links does not contain repeated data
  !

  do iexp = 1, nexp
    do i = 1, experiment(iexp)%nobs - 1
      do j = i + 1, experiment(iexp)%nobs
        if ( experiment(iexp)%observed(i) .eq. experiment(iexp)%observed(j) ) then
          write(*,*) ' ERROR: An observed link is listed repeatedly in experiment: ', &
                     trim(experiment(iexp)%name)
          write(*,*) ' Repeated observed link: ', trim(print_obs(experiment(iexp)%observed(i)))
          stop
        end if
      end do
    end do
  end do

  !
  ! Read atom information from PDB file
  ! 

  open(10,file=pdbfile,status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open PDB file: ', trim(adjustl(pdbfile))
    stop
  end if
  natoms = 0
  do
    read(10,"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    call strclean(record)
    if ( keyword(record) == endread ) exit
    readatom=read_atom(record,error)
    if ( error ) cycle
    if ( .not. isprotein(readatom) ) cycle
    if ( readatoms == 1 ) then
      if ( ishydrogen(readatom) ) cycle
    else if ( readatoms == 2 ) then
      if ( readatom%name /= 'N' .and. &
           readatom%name /= 'CA' .and. &
           readatom%name /= 'C' .and. &
           readatom%name /= 'O' ) cycle
    else if ( readatoms == 3 ) then 
        if ( readatom%name /= 'N' .and. &
             readatom%name /= 'CA' .and. &
             readatom%name /= 'C' .and. &
             readatom%name /= 'CB' .and. &
             readatom%name /= 'O' ) cycle
    end if
    natoms = natoms + 1
  end do
  rewind(10)
  write(*,*)
  write(*,intout) ' Number of atoms read from PDB file: ', natoms
  allocate(atom(natoms),skip(natoms))

  i = 0
  j = 0
  nchains = 1
  do
    read(10,"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    call strclean(record)
    if ( keyword(record) == endread ) exit
    readatom=read_atom(record,error)
    if ( error ) cycle
    if ( .not. isprotein(readatom) ) cycle
    if ( readatoms == 1 ) then
      if ( ishydrogen(readatom) ) cycle
    else if ( readatoms == 2 ) then
      if ( readatom%name /= 'N' .and. &
           readatom%name /= 'CA' .and. &
           readatom%name /= 'C' .and. &
           readatom%name /= 'O' ) cycle
    else if ( readatoms == 3 ) then 
        if ( readatom%name /= 'N' .and. &
             readatom%name /= 'CA' .and. &
             readatom%name /= 'C' .and. &
             readatom%name /= 'CB' .and. &
             readatom%name /= 'O' ) cycle
    end if
    i = i + 1
    atom(i) = readatom
    if ( trim(atom(i)%residue%chain) == "0" ) then
      if ( j == 0 ) then
        j = 1
        write(*,*) " WARNING: Some atoms do not have a chain specifier. Will attribute chain = '0' "
        warning = .true.
      end if
    end if
    atom(i)%index = i
    if ( i > 1 ) then
      if ( atom(i)%residue%chain /= atom(i-1)%residue%chain ) nchains = nchains + 1
    end if
  end do
  close(10)

  ! If interchain and there is only one chain, report error and stop

  write(*,*) ' Number of chains found in PDB file: ', nchains
  if ( interchain ) then
    if ( nchains == 1 ) then
      write(*,*) ' ERROR: Only one chain was found in PDB file, and the interchain keyword was used. '
      stop
    end if
  end if

  ! For each residue of each atom, finds the first at last atoms

  k = -1
  i1 = 1
  do i = 1, natoms
    ii = atom(i)%residue%index
    if ( ii == k ) then
      atom(i)%residue%firstatom = i1
    else
      k = ii
      i1 = i
      atom(i)%residue%firstatom = i1
    end if
  end do

  ! For each atom, finds the last atom of the residue it belongs

  k = -1
  i1 = natoms
  do i = natoms, 1, -1
    ii = atom(i)%residue%index
    if ( ii == k ) then
      atom(i)%residue%lastatom = i1
    else
      k = ii
      i1 = i
      atom(i)%residue%lastatom = i1
    end if
  end do
  
  !
  ! Checking the validity of the input data
  !

  ! Checking if a linktype was defined with ambiguous reactivity 

  error = .false.
  do iexp = 1, nexp
    do i = 1, experiment(iexp)%ntypes
      if ( experiment(iexp)%linktype(i)%atom1%residue%name == &
           experiment(iexp)%linktype(i)%atom2%residue%name ) then
        if ( experiment(iexp)%linktype(i)%atom1%name /= &
             experiment(iexp)%linktype(i)%atom2%name ) then
          write(*,*) ' ERROR: Two different atoms of the same residue define a linktype. '
          write(*,*) '        The reactivity becomes ambiguous. '
          write(*,*) '        Experiment: ', trim(experiment(iexp)%name)
          write(*,*) '        Link type: ', trim(print_linktype(experiment(iexp)%linktype(i)))
          stop
        end if
      end if
    end do
  end do
  if ( error ) stop

  error = .false.
  do iexp = 1, nexp

    checkobs1 : do i = 1, experiment(iexp)%nobs
      do j = 1, natoms
        if ( atom(j) .in. experiment(iexp)%observed(i)%residue1 ) cycle checkobs1
      end do
      write(*,*)
      write(*,*) ' ERROR: First residue of observed link could not be found on the structure. '
      write(*,*) '        Experiment: ', trim(experiment(iexp)%name)
      write(*,*) '        Observed link: ', trim(print_obs(experiment(iexp)%observed(i)))
      error = .true.
    end do checkobs1

    checkobs2 : do i = 1, experiment(iexp)%nobs
      do j = 1, natoms
        if ( atom(j) .in. experiment(iexp)%observed(i)%residue2 ) cycle checkobs2 
      end do
      write(*,*)
      write(*,*) ' ERROR: Second residue of observed link could not be found on the structure. '
      write(*,*) '        Experiment: ', trim(experiment(iexp)%name)
      write(*,*) '        Observed link: ', trim(print_obs(experiment(iexp)%observed(i)))
      error = .true.
    end do checkobs2

    checkdeadend : do i = 1, experiment(iexp)%ndeadends
      do j = 1, natoms
        if ( atom(j) .in. experiment(iexp)%deadend(i) ) cycle checkdeadend
      end do
      write(*,*)
      write(*,*) ' ERROR: Residue of observed deadend could not be found on the structure. '
      write(*,*) '        Experiment: ', trim(experiment(iexp)%name)
      write(*,*) '        Observed deadend: ', trim(print_deadend(experiment(iexp)%deadend(i)))
      error = .true.
    end do checkdeadend

    checktypes1 : do i = 1, experiment(iexp)%ntypes
      do j = 1, natoms
        if ( atom(j) .in. experiment(iexp)%linktype(i)%atom1%residue ) cycle checktypes1
      end do
      write(*,*)
      write(*,*) ' WARNING: First atom of link type does not correspond to any atom of the structure. '
      write(*,*) '          Experiment: ', trim(experiment(iexp)%name)
      write(*,*) '          Link type: ', trim(print_linktype(experiment(iexp)%linktype(i)))
      warning = .true.
    end do checktypes1
 
    checktypes2 : do i = 1, experiment(iexp)%ntypes
      do j = 1, natoms
        if ( atom(j) .in. experiment(iexp)%linktype(i)%atom2%residue ) cycle checktypes2
      end do
      write(*,*)
      write(*,*) ' WARNING: Second atom of link type does not correspond to any atom of the structure. '
      write(*,*) '          Experiment: ', trim(experiment(iexp)%name)
      write(*,*) '          Link type: ', trim(print_linktype(experiment(iexp)%linktype(i)))
      warning = .true.
    end do checktypes2

    checkobstype : do i = 1, experiment(iexp)%nobs
      do j = 1, experiment(iexp)%ntypes
        if ( experiment(iexp)%observed(i) .matches. experiment(iexp)%linktype(j) ) then
          experiment(iexp)%observed(i)%type = j
          cycle checkobstype
        end if
      end do
      write(*,*)
      write(*,*) ' ERROR: Observed link does not match any link type. '
      write(*,*) '        Experiment: ', trim(experiment(iexp)%name)
      write(*,*) '        Obseved link: ', trim(print_obs(experiment(iexp)%observed(i)))
      error = .true.
    end do checkobstype

  end do
  if ( error ) stop

  ! Checking whether atoms that define observed links are missing in the structure

  error = .false.
  do iexp = 1, nexp
    do k = 1, experiment(iexp)%nobs

      ! Type of linker of this observation

      j = experiment(iexp)%observed(k)%type

      ! Checking first residue of observed link (the next question is ok because the 
      ! observed link was already tested for consistency with some linktype previously)

      testatom = experiment(iexp)%linktype(j)%atom1
      if ( experiment(iexp)%observed(k)%residue1%name /= testatom%residue%name ) then
        testatom = experiment(iexp)%linktype(j)%atom2
      end if
      error = .true.
      do i = 1, natoms
        if ( atom(i)%residue%chain == experiment(iexp)%observed(k)%residue1%chain ) then
          if ( atom(i)%residue%index == experiment(iexp)%observed(k)%residue1%index ) then
            if ( atom(i)%residue%name == experiment(iexp)%observed(k)%residue1%name ) then
              if ( atom(i)%name == testatom%name ) then
                error = .false.
                exit
              end if
            end if
          end if
        end if
      end do
      if ( error ) then
        write(*,*) ' ERROR: Atom missing in the structure is required for observed link: '
        write(*,*) '        Missing atom: ', testatom%residue%name, &
                   experiment(iexp)%observed(k)%residue1%chain, &
                   experiment(iexp)%observed(k)%residue1%index, &
                   testatom%name
        stop
      end if

      ! Checking second residue of observed link

      testatom = experiment(iexp)%linktype(j)%atom1
      if ( experiment(iexp)%observed(k)%residue2%name /= testatom%residue%name ) then
        testatom = experiment(iexp)%linktype(j)%atom2
      end if
      error = .true.
      do i = 1, natoms
        if ( atom(i)%residue%chain == experiment(iexp)%observed(k)%residue2%chain ) then
          if ( atom(i)%residue%index == experiment(iexp)%observed(k)%residue2%index ) then
            if ( atom(i)%residue%name == experiment(iexp)%observed(k)%residue2%name ) then
              if ( atom(i)%name == testatom%name ) then
                error = .false.
                exit
              end if
            end if
          end if
        end if
      end do
      if ( error ) then
        write(*,*) ' ERROR: Atom missing in the structure is required for observed link: '
        write(*,*) '        Missing atom: ', testatom%residue%name, &
                   experiment(iexp)%observed(k)%residue2%index, &
                   testatom%name
        stop
      end if
    end do
  end do
 
  ! Writing the observed and deadend data back

  do iexp = 1, nexp

    write(*,*)
    write(*,dashes)
    write(*,*) ' Experiment: ', trim(experiment(iexp)%name)

    ! Printing observed links for this experiment
       
    if ( experiment(iexp)%nobs > 0 ) then
      write(*,*) ' Observed links: '
      do i = 1, experiment(iexp)%nobs
        write(*,"( i5, a )") i, trim(print_obs(experiment(iexp)%observed(i)))
      end do
    end if

   ! Printing the observed deadend data
      
    if ( experiment(iexp)%ndeadends > 0 ) then
      write(*,*) ' Observed deadends: '
      do i = 1, experiment(iexp)%ndeadends
        write(*,*) trim(print_deadend(experiment(iexp)%deadend(i)))
      end do
    end if

  end do
  write(*,*)
  write(*,dashes) 
  write(*,dashes) 

  ! Counting how many atoms in the structure are "reactive", according to link types, save
  ! their indices in reactiveatom vector

  i1 = 0
  reactive1 : do i = 1, natoms
    do iexp = 1, nexp
      do j = 1, experiment(iexp)%ntypes
        if ( atom(i) .in. experiment(iexp)%linktype(j) ) then
          i1 = i1 + 1
          cycle reactive1
        end if
      end do
    end do
  end do reactive1
  write(*,*)
  write(*,intout) ' Number of atoms that might be involved in reactions, according link types: ', i1
  natreactive = i1
  allocate( reactiveatom(natreactive) )
  i1 = 0
  reactive2 : do i = 1, natoms
    do iexp = 1, nexp
      do j = 1, experiment(iexp)%ntypes
        if ( atom(i) .in. experiment(iexp)%linktype(j) ) then
          i1 = i1 + 1
          reactiveatom(i1) = i
          cycle reactive2
        end if
      end do
    end do
  end do reactive2

  ! Counting how many atoms in the structure are "reactive", according to observed reactivity

  write(*,*)
  write(*,*) ' All reactive atoms, according to observations: ' 
  i2 = 0
  reactive3 : do i = 1, natreactive
    i1 = reactiveatom(i)
    do iexp = 1, nexp
      do j = 1, experiment(iexp)%nobs
        if ( atom(i1) .in. experiment(iexp)%observed(j) ) then
          i2 = i2 + 1
          write(*,*) i2, print_atom(atom(i1))
          cycle reactive3
        end if
      end do
      do j = 1, experiment(iexp)%ndeadends
        if ( atom(i1) .in. experiment(iexp)%deadend(j) ) then
          i2 = i2 + 1
          write(*,*) i2, print_atom(atom(i1))
          cycle reactive3
        end if
      end do
    end do
  end do reactive3

  ! Maximum number of pairs of reactive atoms, according to link type definitions

  npairs = natreactive*(natreactive-1)/2
  allocate( observed(npairs,nexp), iobserved(npairs,nexp), &
            obs_reactive(npairs,nexp), type_reactive(npairs,nexp), &
            dlink(npairs,nexp), link(npairs) ) 
  do i = 1, npairs
    allocate(link(i)%exp(nexp))
    link(i)%n_type_expected = 0
    link(i)%n_obs_expected = 0
    link(i)%n_type_consistent = 0
    link(i)%n_obs_consistent = 0
    do iexp = 1, nexp
      link(i)%exp(iexp)%observed = .false.
      link(i)%exp(iexp)%type_reactive = .false.
      link(i)%exp(iexp)%obs_reactive = .false.
      link(i)%exp(iexp)%type_consistent = .false.
      link(i)%exp(iexp)%obs_consistent = .false.
    end do
  end do

  ! Annotating which pairs are of each type, for each experiment

  do iexp = 1, nexp

    i1 = 1
    i2 = 1
    do ii = 1, npairs

      if ( i2 == natreactive ) then
        i1 = i1 + 1
        i2 = i1
      end if
      i2 = i2 + 1

      i = reactiveatom(i1)
      j = reactiveatom(i2)

      link(ii)%atom1 = atom(i)
      link(ii)%atom2 = atom(j)
       
      observed(ii,iexp) = .false.
      obs_reactive(ii,iexp) = .false.
      type_reactive(ii,iexp) = .false.
      dlink(ii,iexp) = 0.d0

      ! Annotating if this pair was observed in this experiment

      do k = 1, experiment(iexp)%nobs
        if ( link(ii) .matches. experiment(iexp)%observed(k) ) then
          observed(ii,iexp) = .true.
          iobserved(ii,iexp) = k
          link(ii)%exp(iexp)%observed = .true.
        end if
      end do
      
      ! Annotating if this pair matches a link type of this experiment

      do k = 1, experiment(iexp)%ntypes
        if ( link(ii) .matches. experiment(iexp)%linktype(k) ) then
          if( type_reactive(ii,iexp) ) then
            write(*,*) ' ERROR: A reactive pair of atoms corresponds to two different types of links '
            write(*,*) '        of the same experiment. The link types are not correctly defined. '
            write(*,*) ' Experiment: ', trim(experiment(iexp)%name)
            write(*,*) ' Atoms: ', print_atom(link(ii)%atom1), print_atom(link(ii)%atom2)
            stop
          else
            type_reactive(ii,iexp) = .true.
            dlink(ii,iexp) = experiment(iexp)%linktype(k)%dist
            link(ii)%exp(iexp)%type_reactive = .true.
            link(ii)%n_type_expected = link(ii)%n_type_expected + 1
          end if
        end if
      end do

      ! Annotating the reactive atom pairs according to this experiment

      r1 = .false.
      r2 = .false.
      do k = 1, experiment(iexp)%nobs
        if ( atom(i) .in. experiment(iexp)%observed(k) ) r1 = .true.
        if ( atom(j) .in. experiment(iexp)%observed(k) ) r2 = .true.
      end do
      do k = 1, experiment(iexp)%ndeadends
        if ( atom(i) .in. experiment(iexp)%deadend(k) ) r1 = .true.
        if ( atom(j) .in. experiment(iexp)%deadend(k) ) r2 = .true.
      end do
      if ( r1 .and. r2 ) then
        linktest%atom1 = atom(i)
        linktest%atom2 = atom(j)
        do k = 1, experiment(iexp)%ntypes
          if ( linktest .matches. experiment(iexp)%linktype(k) ) then 
            obs_reactive(ii,iexp) = .true.
            link(ii)%exp(iexp)%obs_reactive = .true.
            link(ii)%n_obs_expected = link(ii)%n_obs_expected + 1
          end if
        end do
      end if

    end do

  end do

  ! Compute the number of reactive pairs of each experiment, according to the link types
  ! and sequence, according to the observed reactivity

  do iexp = 1, nexp
    j = 0
    experiment(iexp)%nreactive_type = 0
    experiment(iexp)%nreactive_obs = 0 
    do i = 1, npairs
      if ( type_reactive(i,iexp) ) then
        j = j + 1
        experiment(iexp)%nreactive_type = experiment(iexp)%nreactive_type + 1
      end if
      if ( obs_reactive(i,iexp) ) then
        experiment(iexp)%nreactive_obs = experiment(iexp)%nreactive_obs + 1
      end if
    end do
  end do
  write(*,*) 
  write(*,*) ' Number of reactive pairs according to sequence: '
  i1 = 0
  do iexp = 1, nexp
    write(*,"( tr3, a, a, a, i5 )") ' Experiment ', trim(experiment(iexp)%name), ': ', &
            experiment(iexp)%nreactive_type
    i1 = i1 + experiment(iexp)%nreactive_type
  end do
  write(*,*) ' Total number of unique reactive pairs, according to link types: ', i1
  write(*,*) 
  write(*,*) ' Number of reactive pairs according to observations: '
  i1 = 0
  do iexp = 1, nexp
    write(*,"( tr3, a, a, a, i5 )") ' Experiment ', trim(experiment(iexp)%name), ': ', &
            experiment(iexp)%nreactive_obs
    i1 = i1 + experiment(iexp)%nreactive_obs
  end do
  write(*,*) ' Total number of unique reactive pairs, according to observations: ', i1

  ! Printing the list of reactive pairs of atoms

  write(*,*)
  write(*,*) ' List of reactive atom pairs, according to observations: ' 
  j = 0
  do i = 1, npairs
    r1 = .false.
    record = " in experiments: "
    do iexp = 1, nexp
      if ( obs_reactive(i,iexp) ) then
        if ( .not. r1 ) then
          j = j + 1
          r1 = .true. 
          record = trim(record)//trim(experiment(iexp)%name)
        else
          record = trim(record)//','//trim(experiment(iexp)%name)
        end if
      end if
    end do
    if ( r1 ) then
       write(*,"( i5, a, a, a )") j, print_atom(link(i)%atom1), &
                                     print_atom(link(i)%atom2), &
                                     trim(record)
    end if
  end do

  !
  ! Now, starting to set computations
  !

  ! Determine the maximum number of variables possible from linker lengths

  nmax = 0
  do i = 1, npairs
    link(i)%dmaxlink = 0.d0
    link(i)%nbeads = 0
    link(i)%observed = .false.
    link(i)%obs_reactive = .false.
    link(i)%type_reactive = .false.
    if ( interchain .and. &
         ( link(i)%atom1%residue%chain == link(i)%atom2%residue%chain ) ) cycle
    do iexp = 1, nexp
      if ( link(i)%exp(iexp)%type_reactive ) then
        link(i)%dmaxlink = max(link(i)%dmaxlink,dlink(i,iexp))
      end if
      if ( link(i)%exp(iexp)%observed ) link(i)%observed = .true.
      if ( link(i)%exp(iexp)%obs_reactive ) link(i)%obs_reactive = .true.
      if ( link(i)%exp(iexp)%type_reactive ) link(i)%type_reactive = .true.
    end do
    link(i)%nbeads = int(link(i)%dmaxlink/dbond)+1
    link(i)%topodist = -1.d0
    nmax = max(nmax,link(i)%nbeads)
  end do
  nmax = 3*nmax
  allocate( x(nmax), g(nmax), xbest(nmax), sigma(nmax) )

  ! Compute maximum and minimum distances of links, according to observations

  do i = 1, npairs
    link(i)%dmin = 0.d0
    link(i)%dmax = link(i)%dmaxlink
    do iexp = 1, nexp
      if ( link(i)%exp(iexp)%observed ) then
        link(i)%dmax = dmin1(link(i)%dmax,dlink(i,iexp))
      end if
    end do
    do iexp = 1, nexp
      if ( link(i)%exp(iexp)%type_reactive .and. &
           .not. link(i)%exp(iexp)%observed ) then
        if ( dlink(i,iexp) <= link(i)%dmax ) then
          link(i)%dmin = dmax1(link(i)%dmin,dlink(i,iexp))
        end if
      end if
    end do
  end do

  ! Initialize vector containing coordinates only, for computations

  allocate( coor(natoms,3) )
  do i = 1, natoms
    coor(i,1) = atom(i)%x
    coor(i,2) = atom(i)%y
    coor(i,3) = atom(i)%z
  end do

  ! Initialize linked cells

  call initcells()
  if ( print > 0 ) then
    write(*,"( t2,a, 6(f8.3) )") ' Minimum and maximum coordinates of the structure: ',&
                                 xmin, ymin, zmin, xmax, ymax, zmax
  end if

  ! Check which residues are solvent accessible

  call solventaccess(atom)
 
  !
  ! Start computation of topological distances
  !

  ! Print title of LINK output if print == 0

  write(*,*)
  if ( print == 0 ) call printdata(-1,link(1))

  allpairs : do i = 1, npairs

    ! Cycle if this topological link must not be computed

    if ( compute == 1 .and. .not. link(i)%observed ) cycle
    if ( compute == 2 .and. .not. link(i)%obs_reactive ) cycle
    if ( .not. link(i)%type_reactive ) cycle

    ! If only interchain links are to be computed
    if ( interchain .and. &
         ( link(i)%atom1%residue%chain == link(i)%atom2%residue%chain ) ) cycle

    atom1 = link(i)%atom1%index
    atom2 = link(i)%atom2%index
    first(1) = link(i)%atom1%residue%firstatom
    first(2) = link(i)%atom2%residue%firstatom
    last(1) = link(i)%atom1%residue%lastatom
    last(2) = link(i)%atom2%residue%lastatom
    nlinkatoms = link(i)%nbeads

    ! Solvent accessibility of atoms and residues

    link(i)%atom1%accessible = atom(atom1)%accessible
    link(i)%atom1%residue%accessible = atom(atom1)%residue%accessible
    link(i)%atom2%accessible = atom(atom2)%accessible
    link(i)%atom2%residue%accessible = atom(atom2)%residue%accessible

    if ( print > 0 ) then
      write(*,"( a, a, 3(tr2,f8.3) )") ' Reference atom 1: ', &
                    print_atom(link(i)%atom1), &
                    link(i)%atom1%x, link(i)%atom1%y, link(i)%atom1%z
      write(*,"( a, a, 3(tr2,f8.3) )") ' Reference atom 2: ', &
                    print_atom(link(i)%atom2), &
                    link(i)%atom2%x, link(i)%atom2%y, link(i)%atom2%z
      if (link(i)%observed ) write(*,*) ' This link was experimentally observed '
      write(*,intout2) ' First and last atoms of residue 1: ', first(1), last(1)
      write(*,intout2) ' First and last atoms of residue 2: ', first(2), last(2)
    end if

    ! Try to read the properties of this link from a previous log file, if required

    link(i)%status = -1
    if ( readlog /= "none" ) then
      do i1 = 1, nloglines
        record = logline(i1)
        if ( record(3:7) /= "LINK:" ) cycle
        linktest = read_link(logline(i1))
        if ( linktest .eq. link(i) ) then
          if ( linktest%status /= -1 ) then
            link(i)%found = linktest%found
            link(i)%euclidean = linktest%euclidean
            link(i)%topodist = linktest%topodist
            ! If the link was not found previously in the structure, it might have
            ! to be searched for again if the previous dmax was smaller than the
            ! new one
            if ( linktest%status == 3 .or. &
                 linktest%status == 4 .or. &
                 linktest%status == 7 .or. &
                 linktest%status == 8 ) then
              if ( linktest%dmax < link(i)%dmaxlink .and. &
                   linktest%euclidean < link(i)%dmaxlink ) then
                link(i)%status = -1
                exit
              end if
            end if
            ! Set link status, according to the link data
            link(i)%status = linkstatus(link(i))
            exit
          end if
        end if
      end do
    end if

    ! Otherwise, compute the link data from the structure

    if ( link(i)%status == -1 ) then

      ! Set found link to false

      link(i)%found = .false.

      ! Compute the topological distance

      link(i)%euclidean =  dsqrt( (coor(atom2,1) - coor(atom1,1))**2 + & 
                                  (coor(atom2,2) - coor(atom1,2))**2 + & 
                                  (coor(atom2,3) - coor(atom1,3))**2 )
      if ( print > 0 ) write(*,floatout) ' Euclidean distance: ', link(i)%euclidean
      if ( link(i)%euclidean > link(i)%dmaxlink ) then
        if ( printnotfound ) then
          link(i)%status = linkstatus(link(i))
          call linkconsistency(link(i),nexp,experiment)
          call printdata(print,link(i))
        end if
        cycle allpairs
      end if

      ! If some of the residues involved are not accessible to solvent, cycle

      if ( .not. link(i)%atom1%residue%accessible .or. &
           .not. link(i)%atom2%residue%accessible ) then
        if ( printnotfound ) then
          link(i)%status = linkstatus(link(i))
          call linkconsistency(link(i),nexp,experiment)
          call printdata(print,link(i))
        end if
        cycle allpairs
      end if

      ! Number of variables of the optimization problem
    
      n = nlinkatoms*3
    
      if ( print > 0 ) then
        write(*,intout2) ' Reference atoms: ', atom1, atom2
        write(*,*) ' Solvent accessibility of these atoms: ', atom(atom1)%accessible, &
                                                              atom(atom2)%accessible
        write(*,*) ' Solvent accessibility of the residues: ', atom(atom1)%residue%accessible, &
                                                               atom(atom2)%residue%accessible
        write(*,intout) ' Number of atoms of the linker: ', nlinkatoms
        write(*,intout) ' Number of variables of the optimization problem: ', n
        write(*,floatout) ' Force constante for link bonds: ', kbond
      end if

      ! Define the atoms that must be invisible in this path calculation

      call setskip(link(i),atom)

      ! Define the vdw radii of the linker atoms, which changes with mimic chain option
      
      call setsigma(link(i),mimicchain)

      ! Start iterative procedure of linker length minimization

      dpath_best = 1.d30
      best_repeat = 0
      itrial = 0
      do while( best_repeat < nbest .and. itrial < ntrial )
        itrial = itrial + 1

        ! Initial link guess

        call initguess(n,x,iguess)  

        ! Test analytical gradient (debugging purposes only)
        ! call test_grad(n,x,g,computef,computeg)
        ! stop

        ! Minimize the energy of the linker
    
        call callcgnewton(n,x,f,0,computef,computeg,optpars)

        ! Remove any possible remaning overlap
    
        kvdw=50.d0*kvdw
        call callcgnewton(n,x,f,0,computef,computeg,optpars)
        kvdw=kvdwini

        dpath = computedpath(n,x)
        overviol = overlap(n,x)/nlinkatoms
        if ( dmin_maxviol < 1.d-1*vdwrad .and. &
             dbond_maxviol < 1.d-2*dbond ) then

          link(i)%found = .true.
          if ( abs(dpath_best-dpath)/dpath_best < 1.d-2 ) then
            best_repeat = best_repeat + 1
          else
            if ( (dpath_best-dpath)/dpath_best > 1.d-2 ) best_repeat = 0
          end if
          if ( dpath <= dpath_best ) then
            dpath_best = dpath
            do j = 1, n
              xbest(j) = x(j)
            end do
          end if
          if ( print > 0 ) then
            write(*,"(' Trial ', i5, ' Valid path with length = ', &
                  &f8.3, '( overlap = ', f12.5, ' dmin_maxviol = ', f8.3,' )', i5)") &
                  itrial, dpath, overviol, dmin_maxviol, best_repeat
          end if
          if ( quitgood .and. dpath <= link(i)%dmaxlink ) exit
        else
          if ( print > 0 ) then
            write(*,"( ' Trial ', i5, ' Invalid path with length = ', &
                       &f8.3, '( overlap = ', f12.5, ' dmin_maxviol = ', f8.3,' )' )")&
                    itrial, dpath, overviol, dmin_maxviol
          end if
        end if
      end do
      if ( link(i)%found ) then
        link(i)%topodist = dpath_best
      end if
      link(i)%status = linkstatus(link(i))

    end if

    ! Report result of link search for this pair
  
    call linkconsistency(link(i),nexp,experiment)

    if ( link(i)%found ) then

      call printdata(print,link(i))

      ! Write the best link obtained to output file
   
      if ( printlinks ) then
        linkfile = pdbfile
        call cleanname(linkfile)
        write(char1,"( i4 )") link(i)%atom1%residue%index
        write(char2,"( i4 )") link(i)%atom2%residue%index
        linkfile=trim(adjustl(linkdir(1:length(linkdir))))//linkfile(1:length(linkfile)-4)//&
          '_'//trim(adjustl(link(i)%atom1%residue%name))//trim(adjustl(link(i)%atom1%residue%chain))//&
          trim(adjustl(char1))//trim(adjustl(link(i)%atom1%name))//&
          '-'//trim(adjustl(link(i)%atom2%residue%name))//trim(adjustl(link(i)%atom2%residue%chain))//&
          trim(adjustl(char2))//trim(adjustl(link(i)%atom2%name))//'.pdb'
        open(10,file=linkfile,iostat=ioerr)
        if ( ioerr /= 0 ) then
          write(*,*) ' ERROR: Could not create link PDB file: ', trim(adjustl(linkfile))
          write(*,*) '        Perhaps the output directory does not exist.'
          write(*,*) '        Output directory: ', trim(adjustl(linkdir))
          stop
        end if
        write(10,"( 'REMARK F: ', 3(tr2,f12.5) )") f, overlap(n,x), stretch(n,x)
        writeatom = atom(atom1)
        writeatom%residue%name = "LINK"
        writeatom%residue%chain = "A"
        writeatom%residue%index = 1
        writeatom%index = 1
        write(10,"(a)") trim(print_pdbhetatm(writeatom))
        do j = 1, nlinkatoms
          ix = (j-1)*3 + 1
          iy = ix + 1
          iz = ix + 2
          writeatom%index = j + 1 
          writeatom%name = "O"
          writeatom%residue%name = "LINK"
          writeatom%residue%chain = "A"
          writeatom%residue%index = 1
          writeatom%x = xbest(ix)
          writeatom%y = xbest(iy)
          writeatom%z = xbest(iz)
          write(10,"(a)") trim(print_pdbhetatm(writeatom))
        end do
        writeatom = atom(atom2)
        writeatom%residue%name = "LINK"
        writeatom%residue%chain = "A"
        writeatom%residue%index = 1
        writeatom%index = nlinkatoms + 2
        write(10,"(a)") trim(print_pdbhetatm(writeatom))
        close(10)
      end if

    else

      ! If a topological distance was not found for this pair
      if ( printnotfound ) call printdata(print,link(i))

    end if

  end do allpairs
  write(*,dashes)
  if ( readlog /= 'none' ) close(20)
 
  write(*,*) 
  write(*,*) ' RESULTS: '
  
  write(*,*)
  write(*,*) ' For each experiment: '

  do iexp = 1, nexp

    experiment(iexp)%ngood = 0
    experiment(iexp)%nbad = 0
    experiment(iexp)%nreach_type = 0
    experiment(iexp)%nreach_obs = 0
    experiment(iexp)%nmiss_type = 0
    experiment(iexp)%nmiss_obs = 0
    experiment(iexp)%score = 0.d0

    do i = 1, npairs

      ! Only reactive pairs, according to link types, of this experiment are of interest

      if ( .not. link(i)%exp(iexp)%type_reactive ) cycle
      if ( compute == 2 .and. .not. link(i)%exp(iexp)%obs_reactive ) cycle

      !
      ! If a topological link was found 
      !

      if ( link(i)%found ) then

        ! If the topological distance is smaller than the link distance ...

        if ( link(i)%topodist <= dlink(i,iexp) ) then
           
          if ( compute == 3 ) then
            experiment(iexp)%nreach_type = experiment(iexp)%nreach_type + 1
          end if
          if ( link(i)%exp(iexp)%obs_reactive ) then
            experiment(iexp)%nreach_obs = experiment(iexp)%nreach_obs + 1
          end if

          ! ... and the link was observed, this link is a good link 

          if ( link(i)%exp(iexp)%observed ) then
            experiment(iexp)%ngood = experiment(iexp)%ngood + 1
            ii = iobserved(i,iexp)
            experiment(iexp)%score = experiment(iexp)%score + &
                                     experiment(iexp)%observed(ii)%score
          end if

          ! ... and the link was NOT observed, this link is a missing link

          if ( .not. link(i)%exp(iexp)%observed ) then
            if ( compute == 3 ) then
              experiment(iexp)%nmiss_type = experiment(iexp)%nmiss_type + 1
            end if
            if ( link(i)%exp(iexp)%obs_reactive ) then
              experiment(iexp)%nmiss_obs = experiment(iexp)%nmiss_obs + 1
            end if
          end if

        end if

        ! If the topological distance is greater than the link distance ...

        if ( link(i)%topodist > dlink(i,iexp) ) then

          ! ... and the link was observed, this observation is bad
      
          if ( link(i)%exp(iexp)%observed ) then
            experiment(iexp)%nbad = experiment(iexp)%nbad + 1
            ii = iobserved(i,iexp)
            experiment(iexp)%score = experiment(iexp)%score - &
                                     experiment(iexp)%observed(ii)%score
          end if

        end if

      !
      ! If the link was not found (meaning that its distance is greater than dmaxlink)
      !
      
      else if ( .not. link(i)%found ) then

        ! If the link was observed, this observation is bad
      
        if ( link(i)%exp(iexp)%observed ) then
          experiment(iexp)%nbad = experiment(iexp)%nbad + 1
        end if

      end if

    end do

    ! Number of pairs of reactive atoms outside linker reach

    if ( compute == 3 ) then
      experiment(iexp)%noutreach_type = experiment(iexp)%nreactive_type - experiment(iexp)%nreach_type
    end if
    experiment(iexp)%noutreach_obs = experiment(iexp)%nreactive_obs - experiment(iexp)%nreach_obs

    !
    ! Likelyhood of the observed result
    !

    experiment(iexp)%pgood = dble(experiment(iexp)%nreach_obs)/experiment(iexp)%nreactive_type
    experiment(iexp)%pbad = dble(experiment(iexp)%nreactive_type - experiment(iexp)%nreach_obs)/&
                            experiment(iexp)%nreactive_type
    experiment(iexp)%likelyhood = lnf(experiment(iexp)%nobs) - &
                                  lnf(experiment(iexp)%ngood) - &
                                  lnf(experiment(iexp)%nbad) + &
                                  nlnp(experiment(iexp)%ngood,experiment(iexp)%pgood) + &
                                  nlnp(experiment(iexp)%nbad,experiment(iexp)%pbad)
    experiment(iexp)%likelyhood = 1.d0 - dexp(experiment(iexp)%likelyhood) 
   
    !
    ! User-Likelyhood: same thing, but using pgood and pbad from user input, if provided
    !

    if ( pbad > 0.d0 .and. pgood > 0.d0 ) then

      experiment(iexp)%userlikelyhood = lnf(experiment(iexp)%nobs) - &
                                        lnf(experiment(iexp)%ngood) - &
                                        lnf(experiment(iexp)%nbad) + &
                                        nlnp(experiment(iexp)%ngood,pgood) + &
                                        nlnp(experiment(iexp)%nbad,pbad)
      experiment(iexp)%userlikelyhood = 1.d0 - dexp(experiment(iexp)%userlikelyhood) 

    end if
   
    ! Report results for this experiment

    write(*,*)  
    write(*,dashes)  
    write(*,*)  
    write(*,*) ' Experiment: ', trim(experiment(iexp)%name)
    write(*,*)
    write(*,intout) '   Number of type-reactive pairs of atoms: ', experiment(iexp)%nreactive_type
    if ( compute == 3 ) then
      write(*,intout) '   Number of type-reactive pairs of atoms within linker reach: ', experiment(iexp)%nreach_type
      write(*,intout) '   Number of type-reactive pairs of atoms outside linker reach: ', experiment(iexp)%noutreach_type
      write(*,intout) '   Missing links, according to the structure and type-reactivity: ', experiment(iexp)%nmiss_type
    end if
    write(*,*)
    write(*,intout) '   Number of observed-reactive pairs of atoms: ', experiment(iexp)%nreactive_obs
    write(*,intout) '   Number of observed-reactive pairs of atoms within linker reach: ', experiment(iexp)%nreach_obs
    write(*,intout) '   Number of observed-reactive pairs of atoms outside linker reach: ', experiment(iexp)%noutreach_obs
    write(*,intout) '   Missing links, according to the structure and observed-reactivity: ', experiment(iexp)%nmiss_obs
    if ( experiment(iexp)%nobs > 0 ) then
      write(*,*)
      write(*,intout) '   Number of observed links: ', experiment(iexp)%nobs
      write(*,intout) '   Number of observed links consistent with the structure: ', experiment(iexp)%ngood
      write(*,intout) '   Number of observed links NOT consistent with the structure: ', experiment(iexp)%nbad
      write(*,*)
      write(*,floatout) '   Sum of scores of observed links: ', experiment(iexp)%score
      write(*,*)
      write(*,floatout) '   Sensitivity of the cross-linking reaction: ', dble(experiment(iexp)%ngood)/experiment(iexp)%nreach_obs
      write(*,floatout) '   False-assignment probability: ', dble(experiment(iexp)%nbad)/experiment(iexp)%nreactive_type
      write(*,*)
      write(*,floatout) '   Likelyhood of the experimental result: ', experiment(iexp)%likelyhood
      write(*,floatout) '   Log-likelyhood of the experimental result: ', dlog(experiment(iexp)%likelyhood)
      if ( pbad > 0.d0 .and. pgood > 0.d0 ) then
        write(*,*)
        write(*,floatout) '   Likelyhood using user-defined pbad and pgood: ', experiment(iexp)%userlikelyhood
        write(*,floatout) '   Log-likelyhood using user-defined pbad and pgood: ', dlog(experiment(iexp)%userlikelyhood)
        write(*,"( t4, 3(a,f8.3))") ' Using: pgood = ', pgood, '; pbad = ', pbad
      end if
   end if
    
  end do

  write(*,*)
  write(*,dashes)
  
  !
  ! Final results
  !
  
  ! Computing the total number of observations (only used for later printing)

  nobs = 0
  ngood = 0
  totscore = 0.d0
  do iexp = 1, nexp
    nobs = nobs + experiment(iexp)%nobs
    ngood = ngood + experiment(iexp)%ngood
    totscore = totscore + experiment(iexp)%score
  end do

  ngooddist = 0
  nbaddist = 0
  nmisslinks = 0
  do i = 1, npairs
    if ( compute == 1 .and. .not. link(i)%observed ) cycle
    if ( compute == 2 .and. .not. link(i)%obs_reactive ) cycle
    if ( .not. link(i)%type_reactive ) cycle
    if ( link(i)%status == 0 .or. &
         link(i)%status == 6 .or. &
         link(i)%status == 7 .or. &
         link(i)%status == 8 ) then
      ngooddist = ngooddist + 1
    else if ( link(i)%status == 1 .or. &
              link(i)%status == 2 .or. &
              link(i)%status == 3 .or. &
              link(i)%status == 4 ) then
      nbaddist = nbaddist + 1
    else if ( link(i)%status == 5 ) then
      nmisslinks = nmisslinks + 1
    end if
  end do

  write(*,*)
  write(*,*) ' FINAL RESULTS: '
  write(*,*)
  write(*,"( t3, a, i5, a )") ' RESULT0: ', ngood, ' : Number of observations that are consistent with the structure.' 
  write(*,*)
  write(*,"( t3, a, i5, a )") ' RESULT1: ', ngooddist, ' : Number of topological distances consistent with all observations. ' 
  write(*,"( t3, a, i5, a )") ' RESULT2: ', nbaddist, ' : Number of topological distances NOT consistent with observations.  ' 
  write(*,"( t3, a, i5, a )") ' RESULT3: ', nmisslinks, ' : Number of links with missing observations.  ' 

  if ( nobs > 0 ) then
    likelyhood = 1.d0
    userlikelyhood = 1.d0
    do iexp = 1, nexp
      likelyhood = likelyhood*experiment(iexp)%likelyhood 
      if ( pbad > 0.d0 .and. pgood > 0.d0 ) then
        userlikelyhood = userlikelyhood*experiment(iexp)%userlikelyhood
      end if
    end do
    write(*,*)
    write(*,"( t3, a, f12.5, a )") ' RESULT4: ', totscore, ' : Sum of scores of observed links of all experiments. '
    write(*,*)
    write(*,"( t3, a, f12.5, a )") ' RESULT5: ', likelyhood, ' : Likelyhood of the set of experimental results. '
    write(*,"( t3, a, f12.5, a )") ' RESULT6: ', dlog(likelyhood), ' : Log-likelyhood of the set of experimental results. '
    if ( pbad > 0.d0 .and. pgood > 0.d0 ) then
      write(*,*)
      write(*,"( t3, 3(a,f8.3))") ' Using: pgood = ', pgood, '; pbad = ', pbad
      write(*,"( t3, a, f12.5, a )") ' RESULT7: ', userlikelyhood, ' : Likelyhood of the set of experimental results. '
      write(*,"( t3, a, f12.5, a )") ' RESULT8: ', dlog(userlikelyhood), ' : Log-likelyhood of the set of experimental results. '
    end if
  end if

  if ( warning ) then
    write(*,*)
    write(*,dashes)
    write(*,*) ' ATTENTION: There are WARNINGS of input file parsing. Please check the output carefully. '
    write(*,dashes)
    write(*,*)
  end if
   
  write(*,*)
  write(*,hashes)

end program topolink



