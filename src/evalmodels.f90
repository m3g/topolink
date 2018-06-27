!
! Program EvalModels
!
! Reads a series of TopoLink output files and a file containing the scores
! (for example, a LovoAlign alignment file, or any other score), and outputs
! the statistics of link satisfaction obtained by TopoLink as a function of
! this score. 
!
! This is to verify in any of the statistics of links satisfaction is
! correlated with the quality of the score chosen.
!
! L. Martinez
! Instiute of Chemistry - University of Campinas
! September 9, 2016
! http://leandro.iqm.unicamp.br/topolink
!
program evalmodels

  use ioformat 
  use topolink_data
  use topolink_operations
  use string_operations

  implicit none
  integer :: i, j, ilink, imodel
  integer :: nargs, nmodels, ioerr, nlinks, scorecol, modelcol, sortcol, nminmax
  character(len=max_string_length) :: loglist, scorelist, record, record2, record3, line, name, output
  logical :: error
  type(specific_link) :: linktemp
  type(modeldata), allocatable :: model(:)
  integer :: model_index

  ! Print title

  call title()
  write(*,*) ' EVALMODELS: Associate TopoLink results with model quality scores. '
  write(*,*)
  write(*,dashes)

  ! Read list of log files from the command line

  nargs = iargc()
  if ( nargs /= 5 .and. nargs /= 6 ) then
    write(*,*)
    write(*,*) ' ERROR: Run with: evalmodels loglist.txt scores.dat output.dat -c[int] -m[int] -s[int]'
    write(*,*)
    write(*,*) ' Where: loglist.txt is the file containing a list of TopoLink logs. '
    write(*,*) '        scores.dat is a LovoAlign log file or a list of model names with scores. '
    write(*,*) '        output.dat is the name of the output file to be created. '
    write(*,*)
    write(*,*) ' The scores are, for a LovoAlign log file: -c3: TM-score '
    write(*,*) '                                           -c5: RMSD of all atoms. '
    write(*,*) '                                           -c8: GDT_TS score. '
    write(*,*) '                                           -c9: GDT_HA score. '
    write(*,*) ' If the score list is not a LovoAlign log file, -c[int] indicates the column of the '
    write(*,*) ' list containing the score. '
    write(*,*)
    write(*,*) ' The -m[int] argument indicates the column of scores.dat containing the model name. '
    write(*,*)
    write(*,*) ' The -s[int] (optional): Order output using the data of this column of output file '
    write(*,*)
    write(*,*) ' More details at: http://leandro.iqm.unicamp/topolink '
    write(*,*)
    write(*,hashes)
    stop
  end if
  call getarg(1,loglist)

  ! Read from which column the scores have to be read from the scorelist

  scorecol = 0
  modelcol = 0
  sortcol = 0
  do i = 4, max(5,nargs)
    call getarg(i,record)
    if ( record(1:2) == "-c" ) then
      read(record(3:length(record)),*) scorecol
    end if
    if ( record(1:2) == "-m" ) then
      read(record(3:length(record)),*) modelcol
    end if
    if ( record(1:2) == "-s" ) then
      read(record(3:length(record)),*) sortcol
    end if
  end do
  if ( scorecol == 0 ) then
    write(*,*) ' ERROR: Column of ', trim(adjustl(scorelist)), ' containing scores not defined. Ex: -c2 '
    stop
  end if
  if ( modelcol == 0 ) then
    write(*,*) ' ERROR: Column of ', trim(adjustl(scorelist)), ' containing model names not defined. Ex: -m1 '
    stop
  end if
  call getarg(3,output)

  ! Open score file (might be a lovoalign log file, or simply a list of names and scores)
  
  call getarg(2,scorelist)
  open(10,file=scorelist,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(adjustl(scorelist))
    stop
  end if

  ! Read number of models
  nmodels = 0
  do
    read(10,string_read,iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) (record2, i=1, modelcol)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Incorrect number of columns in file: ', trim(adjustl(scorelist))
      stop
    end if
    read(record,*,iostat=ioerr) (record3, i=1, scorecol)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Incorrect number of columns in file: ', trim(adjustl(scorelist))
      stop
    end if
    nmodels = nmodels + 1
  end do
  write(*,*) ' Number of models found in data file: ', nmodels
  allocate(model(nmodels))

  ! Read model names and scores from score log file
  rewind(10)
  imodel = 0
  do
    read(10,string_read,iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if ( comment(record) ) cycle
    read(record,*,iostat=ioerr) (record2, i=1, modelcol)
    read(record,*,iostat=ioerr) (record3, i=1, scorecol)
    imodel = imodel + 1
    model(imodel)%name = basename(record2)
    read(record3,*,iostat=ioerr) model(imodel)%score
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not read score for model: ', trim(adjustl(record2))
      write(*,*) '        in file: ', trim(adjustl(scorelist))
      stop
    end if
  end do 
  close(10)

  ! Order all models by name
  call sort_by_name(nmodels,model)

  write(*,*) ' Reading TopoLink log files ... '
  open(10,file=loglist,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open log list file: ', trim(adjustl(loglist))
    stop
  end if
  i = 0
  do 
    read(10,string_read,iostat=ioerr) record
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
      write(*,*) ' ERROR: A model was not found in score file: ', trim(adjustl(loglist))
      write(*,*) '        Model: ', trim(adjustl(name))
      stop
    end if
    i = i + 1
    call progress(i,1,nmodels)

    !
    ! Check the number of links reported in this file
    !
    nlinks = 0
    do 
      read(20,string_read,iostat=ioerr) line
      if ( ioerr /= 0 ) exit
      if ( line(3:7) == "LINK:" ) nlinks = nlinks + 1
    end do
    model(imodel)%nlinks = nlinks
    allocate(model(imodel)%link(nlinks),model(imodel)%linkindex(nlinks))

    !
    ! Read model data
    ! 
    rewind(20)
    ilink = 0
    nminmax = 0 
    do 
      read(20,string_read,iostat=ioerr) line
      if ( ioerr /= 0 ) exit
      if ( line(3:7) == "LINK:" ) then
        linktemp = read_link(line)
        ilink = ilink + 1
        model(imodel)%link(ilink) = linktemp
        if ( linktemp%status == 0 ) then
          if ( linktemp%dmin > 0.d0 .and. &
               linktemp%dmin /= linktemp%dmax ) then
            nminmax = nminmax + 1
          end if
        end if
      end if
      if ( line(4:11) == "RESULT0:") read(line(12:max_string_length),*,iostat=ioerr) model(imodel)%nobscons
      if ( ioerr /= 0 ) model(imodel)%nobscons = 0
      if ( line(4:11) == "RESULT1:") read(line(12:max_string_length),*,iostat=ioerr) model(imodel)%ntopcons
      if ( ioerr /= 0 ) model(imodel)%ntopcons = 0
      if ( line(4:11) == "RESULT2:") read(line(12:max_string_length),*,iostat=ioerr) model(imodel)%ntopnot
      if ( ioerr /= 0 ) model(imodel)%ntopnot = 0
      if ( line(4:11) == "RESULT3:") read(line(12:max_string_length),*,iostat=ioerr) model(imodel)%nmiss
      if ( ioerr /= 0 ) model(imodel)%nmiss = 0
      if ( line(4:11) == "RESULT4:") read(line(12:max_string_length),*,iostat=ioerr) model(imodel)%sumscores
      if ( ioerr /= 0 ) model(imodel)%sumscores = 0.
      if ( line(4:11) == "RESULT5:") read(line(12:max_string_length),*,iostat=ioerr) model(imodel)%likeli
      if ( ioerr /= 0 ) model(imodel)%likeli = 0.
      if ( line(4:11) == "RESULT6:") read(line(12:max_string_length),*,iostat=ioerr) model(imodel)%loglikeli
      if ( ioerr /= 0 ) model(imodel)%loglikeli = 0.
      if ( line(4:11) == "RESULT7:") read(line(12:max_string_length),*,iostat=ioerr) model(imodel)%usrlike
      if ( ioerr /= 0 ) model(imodel)%usrlike = 0.
      if ( line(4:11) == "RESULT8:") read(line(12:max_string_length),*,iostat=ioerr) model(imodel)%usrloglike
      if ( ioerr /= 0 ) model(imodel)%usrloglike = 0.
    end do
    model(imodel)%nminmax = nminmax

    close(20)
  end do
  close(10)

  ! Indexing the links

  write(*,*) ' Indexing links ... '
  imodel = 1
  do i = 1, model(imodel)%nlinks
    model(imodel)%linkindex(i) = i
  end do
  call progress(imodel,1,nmodels)
  do imodel = 2, nmodels 
    call progress(imodel,1,nmodels)
    do i = 1, model(imodel)%nlinks
      jdo : do j = 1, model(1)%nlinks 
        if ( model(imodel)%link(i) .eq. model(1)%link(j) ) then
          model(imodel)%linkindex(i) = model(1)%linkindex(j)
          exit jdo
        end if
      end do jdo
    end do 
  end do

  ! Sorting models according to the user desire

  if ( sortcol > 0 ) then
    write(*,*) ' Sorting models using column output: ', sortcol
    call sort_by_value(nmodels,model,sortcol)
  end if

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
  write(10,"(a)") "# EvalModels output file. " 
  write(10,"(a)") "#"
  write(10,"(a,a)") "# Log file list: ", trim(adjustl(loglist))
  write(10,"(a,a)") "# Score (possibly LovoAlign log) file: ", trim(adjustl(scorelist))
  write(10,"(a,i8)") "# Number of models ", nmodels
  write(10,"(a)") "#"
  write(10,"(a,i5,a)") "# Score: Model quality score, obtained from column ", scorecol,&
                       " of the score file. " 
  write(10,"(a)") "#"
  write(10,"(a)") "# RESULT0: Number of consistent observations. "
  write(10,"(a)") "# RESULT1: Number of topological distances consistent with all observations. "
  write(10,"(a)") "# RESULT2: Number of topological distances NOT consistent with observations. "
  write(10,"(a)") "# RESULT3: Number of missing links in observations. "
  write(10,"(a)") "# RESULT4: Number of distances with min and max bounds that are consistent."
  write(10,"(a)") "# RESULT5: Sum of the scores of observed links in all observations. "
  write(10,"(a)") "# RESULT6: Likelihood of the structural model, based on observations. "
  write(10,"(a)") "#"
  write(10,"(a)") "# More details at: http://leandro.iqm.unicamp.br/topolink"
  write(10,"(a)") "#"
  write(10,"(a)") "#      Score   RESULT0   RESULT1   RESULT2   RESULT3   RESULT4       RESULT5       RESULT6  MODEL"
  do imodel = 1, nmodels
    call progress(imodel,1,nmodels)
    write(10,"( f12.5,5(tr2,i8),tr2,f12.5,tr2,e12.5,tr2,a )") &
                model(imodel)%score, &
                model(imodel)%nobscons, &
                model(imodel)%ntopcons, &
                model(imodel)%ntopnot, &
                model(imodel)%nmiss, &
                model(imodel)%nminmax, &
                model(imodel)%sumscores, &
                model(imodel)%likeli,&
                trim(adjustl(model(imodel)%name))
  end do

  close(10)

  write(*,*) ' Wrote output file: ', trim(adjustl(output))
  write(*,*)
  write(*,hashes)

end program evalmodels









