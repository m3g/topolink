!
! Program LinkEnsemble
!
! Reads the output of G-score for a set of models, and reads
! the output of TopoLink for each model. Compute, for a given
! predefined relative "probability" of the model, that is,
! for relative G-scores, the set of crosslinks
! that are satisfied. This is done by running from higher to
! lower G-score, accumulating the set of links that are 
! satisfied.
!
! L. Martinez
! Instiute of Chemistry - University of Campinas
! September 17, 2016
! http://leandro.iqm.unicamp.br/topolink
!

program linkensemble

  use ioformat 
  use topolink_data
  use topolink_operations
  use string_operations

  implicit none
  integer :: i, j, ilink, imodel, nobserved
  integer :: nargs, nmodels, ioerr, nlinks, nsatisfied
  integer, allocatable :: satisfied(:)
  double precision :: gscore, degree
  character(len=200) :: loglist, gscorefile, record, name, output, line
  logical :: error
  type(specific_link) :: linktemp
  type(modeldata), allocatable :: model(:)
  integer :: model_index

  ! Print title

  call title()
  write(*,*) ' LINKENSEMBLE: Compute the links that are satisfied by the most probable set of models. '
  write(*,*)
  write(*,dashes)

  ! Read list of log files from the command line

  nargs = iargc()
  if ( nargs /= 3 ) then
    write(*,*)
    write(*,*) ' ERROR: Run with: linkensemble loglist.txt gscores.dat output.dat'
    write(*,*)
    write(*,*) ' Where: loglist.txt is the file containing a list of TopoLink logs. '
    write(*,*) '        gscores.dat is a log file of the G-score program.'
    write(*,*) '        output.dat is the name of the output file to be created. '
    write(*,*)
    write(*,*) ' More details at: http://leandro.iqm.unicamp/topolink '
    write(*,*)
    write(*,hashes)
    stop
  end if

  call getarg(1,loglist)
  call getarg(2,gscorefile)
  call getarg(3,output)

  ! Open score file (might be a lovoalign log file, or simply a list of names and scores)
  
  open(10,file=gscorefile,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(gscorefile))
    stop
  end if

  ! Read number of models
  nmodels = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) gscore, degree, name
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Error reading gscore and name from: ', trim(adjustl(gscorefile))
      stop
    end if
    nmodels = nmodels + 1
  end do
  write(*,*) ' Number of models found in gscore file: ', nmodels
  allocate(model(nmodels))

  ! Read model names and scores from score log file
  rewind(10)
  imodel = 0
  do
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) gscore, degree, name
    imodel = imodel + 1
    model(imodel)%name = basename(name)
    model(imodel)%score = gscore
    model(imodel)%degree = degree
  end do 
  close(10)

  ! Sort models by name
  call sort_by_name(nmodels,model)

  write(*,*) ' Reading TopoLink log files ... '
  open(10,file=loglist,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open log list file: ', trim(adjustl(loglist))
    stop
  end if
  i = 0
  do 
    read(10,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    open(20,file=record,status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not open TopoLink log file: ', trim(adjustl(record))
      write(*,*) '        Listed in: ', trim(adjustl(loglist))
      stop
    end if
    !
    ! Check which model is this
    !
    name = basename(record)
    imodel = model_index(name,model,nmodels,error)
    if ( error ) then
    !  write(*,*) ' Warning: A model was not found in list: ', trim(adjustl(loglist))
    !  write(*,*) '          Model: ', trim(adjustl(name))
      cycle
    end if
    i = i + 1
    call progress(i,1,nmodels)

    !
    ! Check the number of links reported in this file
    !

    nlinks = 0
    do 
      read(20,"(a200)",iostat=ioerr) line
      if ( ioerr /= 0 ) exit
      if ( line(3:7) == "LINK:" ) nlinks = nlinks + 1
    end do
    if ( nlinks == 0 ) then
      write(*,*) ' ERROR: No LINK line reported for model: ', model(imodel)%name
      stop
    end if
    model(imodel)%nlinks = nlinks
    allocate(model(imodel)%link(nlinks),model(imodel)%linkindex(nlinks))

    !
    ! Read model data
    ! 

    rewind(20)
    ilink = 0
    do 
      read(20,"(a200)",iostat=ioerr) line
      if ( ioerr /= 0 ) exit
      if ( line(3:7) == "LINK:" ) then
        linktemp = read_link(line)
        ilink = ilink + 1
        model(imodel)%link(ilink) = linktemp
      end if
      if ( line(4:11) == "RESULT0:") read(line(12:200),*) model(imodel)%nobscons
      if ( line(4:11) == "RESULT1:") read(line(12:200),*) model(imodel)%ntopcons
      if ( line(4:11) == "RESULT2:") read(line(12:200),*) model(imodel)%ntopnot
      if ( line(4:11) == "RESULT3:") read(line(12:200),*) model(imodel)%nmiss
      if ( line(4:11) == "RESULT4:") read(line(12:200),*) model(imodel)%sumscores
      if ( line(4:11) == "RESULT5:") read(line(12:200),*) model(imodel)%likely
      if ( line(4:11) == "RESULT6:") read(line(12:200),*) model(imodel)%loglikely
      if ( line(4:11) == "RESULT7:") read(line(12:200),*) model(imodel)%usrlike
      if ( line(4:11) == "RESULT8:") read(line(12:200),*) model(imodel)%usrloglike
    end do

    close(20)
  end do
  close(10)
  if ( i == 0 ) then  
    write(*,*) ' ERROR: No corresponding data was found between files. '
    write(*,*) '        Are you sure the correct files were used? '
    stop
  end if

  ! Order all models by gscore, maybe they are not already sorted

  write(*,*) ' Sorting models by G-score value ... '
  call sort_by_value(nmodels,model,9)

  ! Indexing the links

  imodel = 1
  do i = 1, model(imodel)%nlinks
    model(imodel)%linkindex(i) = i
  end do
  do imodel = 2, nmodels 
    do i = 1, model(imodel)%nlinks
      jdo : do j = 1, model(1)%nlinks 
        if ( model(imodel)%link(i) .eq. model(1)%link(j) ) then
          model(imodel)%linkindex(i) = model(1)%linkindex(j)
          exit jdo
        end if
      end do jdo
    end do 
  end do

  ! Checking the number of observed links

  imodel = 1
  nobserved = 0
  do i = 1, model(imodel)%nlinks
    if ( model(imodel)%link(i)%observed ) then
      nobserved = nobserved + 1
    end if
  end do
  write(*,"(a,i8)") "  Number of observed crosslinks: ", nobserved
  if ( nobserved == 0 ) then
    write(*,*) " ERROR: Cannot run if number of observed crosslinks is zero. "
    stop
  end if
  allocate( satisfied(nobserved) )
  do i = 1, nobserved
    satisfied(i) = 0
  end do

  !
  ! Write output file
  ! 

  write(*,*) ' Writing output file ... '
  open(10,file=output,iostat=ioerr)
  if ( ioerr /= 0 ) then 
    write(*,*) ' ERROR: Could not open output file: ', trim(adjustl(output))
    stop
  end if
  
  write(10,"(a)") "# TopoLink" 
  write(10,"(a)") "#"
  write(10,"(a)") "# LinkEnsemble output file. " 
  write(10,"(a)") "#"
  write(10,"(a,a)") "# Log file list: ", trim(adjustl(loglist))
  write(10,"(a,a)") "# G-score file: ", trim(adjustl(gscorefile))
  write(10,"(a,i8)") "# Number of models ", nmodels
  write(10,"(a,i8)") "# Number of observed crosslinks: ", nobserved
  write(10,"(a)") "#"
  imodel = 1
  do i = 1, model(imodel)%nlinks
    if ( model(imodel)%link(i)%observed ) then
      write(10,"(a,i3,tr1,a)") "# ", i, trim(adjustl(print_link(model(imodel)%link(i))))
    end if
  end do
  write(10,"(a)") "#"
  write(10,"(a)") "# Nmodel: Number of crosslinks satisfied by this model. "
  write(10,"(a)") "# RelatP: Relative probability of this model (G-score ratio to best model)."
  write(10,"(a)") "# DeltaG: RelatP converted to DeltaG (kcal/mol)."
  write(10,"(a)") "# Ntot: Total number of links satisfied by the ensemble up to this model."
  write(10,"(a)") "# Next: link indexes according to list above."
  write(10,"(a)") "#"
  write(record,*) "(a,",nobserved,"(tr1,i3))"
  write(10,record) "#             Model  Nmodel     RelatP       DeltaG  Ntot",(i,i=1,nobserved)
  nsatisfied = 0
  write(record,*) "(i8,tr1,a,tr1,i5,2(tr1,f12.5),tr1,i5,",nobserved,"(tr1,i3))"
  do imodel = 1, nmodels
    model(imodel)%nobsgood = 0
    do i = 1, model(imodel)%nlinks
      ilink = model(imodel)%linkindex(i)
      if ( model(imodel)%link(ilink)%observed ) then
        if ( model(imodel)%link(ilink)%status == 0 ) then
          model(imodel)%nobsgood = model(imodel)%nobsgood + 1
          if ( satisfied(ilink) == 0 ) then
            nsatisfied = nsatisfied + 1
            satisfied(ilink) = 1
          end if
        end if
      end if
    end do
    write(10,record) &
                imodel, &
                trim(adjustl(model(imodel)%name)),&
                model(imodel)%nobsgood,&
                model(imodel)%degree / model(1)%degree, & 
                model(imodel)%score - model(1)%score, & 
                nsatisfied, &
                (satisfied(j),j=1,nobserved)
  end do
 
  write(*,*) ' Wrote output file: ', trim(adjustl(output))
  write(*,*)
  write(*,hashes)

end program linkensemble









