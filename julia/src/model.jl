#
# Function to read the model data from a TopoLink log file and, optionally,
# from a compactlog alignment file
#

function model( logfilename :: String; compactlog :: CompactLog = CompactLog(), 
                                       cutoff = 0.5,
                                       fdist = inverse_minus_one )
                          
  # Read TopoLink log data for this model

  logdata = readlog( logfilename )

  # Read compactlog data for this model, to compute consensus scores

  if compactlog.nmodels != 0

    degree = TopoLink.degree( compactlog, logdata.name, cutoff = cutoff )
    davis = TopoLink.davis( compactlog, logdata.name, f = fdist )
    gscore = -0.593*log( degree / ( compactlog.nmodels - 1. ) )

  else

    gscore = 0.
    degree = 0
    davis = 0.

  end

  return Model( logdata.name, logdata.pdb, logdata.nlinks, logfilename, logdata.link, 
                logdata.nconsist, logdata.nnotcons, logdata.nmissing, gscore, degree, davis )

end

function Base.show( io :: IO, model :: Model )
  if model.degree > 0 
    print( model.name," data with ", model.nlinks," links and consensus score data." )
  else
    print( model.name," data with ", model.nlinks," links without consensus score data." )
  end
end

#
# Functions to read the data of a series of models from a log list, given as a filename
# or an actual list
#

using ProgressMeter

# Read from list of log files (a vector)

function models( loglist :: Vector{String}; compactlog :: CompactLog = CompactLog(), 
                                            cutoff = 0.5,
                                            fdist = inverse_minus_one )

  nlogs = length(loglist)
  models = Vector{Model}(undef,nlogs)

  
  p = Progress(nlogs,1," Reading list of log files: ")
  ilog = 0
  for logfilename in loglist
    ilog = ilog + 1
    models[ilog] = model(logfilename, compactlog=compactlog, cutoff=cutoff, fdist=fdist) 
    update!(p,ilog)
  end

  return models

end

# Read from list of log files given by the name of a file that contains the list

function models( loglistname :: String; compactlog :: CompactLog = CompactLog(), 
                                        cutoff = 0.5,
                                        fdist = inverse_minus_one )

  loglist = open( logliting Pame, "r" )
  nlogs = 0
  for filename in eachline(loglist)
    nlogs = nlogs + 1
  end
  seekstart(loglist)
  list = Vector{String}(undef,nlogs)
  ilog = 0
  for filename in eachline(loglist)
    ilog = ilog + 1
    list[ilog] = filename
  end

  return models( list, compactlog=compactlog, cutoff=cutoff, fdist=fdist )

end

function Base.show( io :: IO, ::MIME"text/plain", models :: Array{Model} )
  print(" Log list with data for ", length(models)," models. ")
end
                          





