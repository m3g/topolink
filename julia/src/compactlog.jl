#
# Functions to read and use the compactlog file
#

using ProgressMeter

#
# Structure that contains the compactlog data
#

struct CompactLog

  file :: String
  nmodels :: Int64
  name :: Vector{String}
  score :: Vector{Float64}

  alignment_file :: String
  pdb_list :: String
  score_type :: String

end

#
# Function to read the compact log file
#

function compactlog( filename :: String )

  local nmodels :: Int64
  local name :: Vector{String}
  local score :: Vector{Float64}
  local alignment_file :: String
  local pdb_list :: String
  local score_type :: String

  file = open(filename,"r")

  seekend(file) ; fileSize = position(file) ; seekstart(file)
  p = Progress(fileSize,1," Reading compactlog file: ")

  imodel = 0
  nmodels = 0
  ipair = 0
  i = 0
  for line in eachline(file) 
    update!(p,position(file))
    rec = split(line)
    if rec[1] == "#"
      if rec[2] == "Alignment"
        alignment_file = strip(line[32:length(line)])
      end
      if rec[2] == "PDB"
        pdb_list = strip(line[12:length(line)])
      end 
      if rec[2] == "Score"
        score_type = strip(line[14:length(line)])
      end
    elseif nmodels == 0 
      # Read number of models and allocate vectors
      nmodels = parse(Int64,line)
      name = Vector{String}(undef,nmodels)
      score = Vector{Float64}(undef,div((nmodels*(nmodels-1)),2))
    elseif nmodels > 0 && imodel < nmodels
      # Read model names
      imodel = imodel + 1
      name[imodel] = rec[2]
    else
      # Read score matrix (upper diagonal, as a vector)
      for data in rec
        i = i + 1
        score[i] = parse(Float64,data)
      end
    end
  end
  close(file)
  
  compactlog = CompactLog(filename,nmodels,name,score,alignment_file,pdb_list,score_type)

  return compactlog
end

CompactLog() = CompactLog("none",0,[],[],"none","none","none")

function Base.show( io :: IO, c :: CompactLog )
  print( " Compact log file with data for ", c.nmodels, " models, score type: ", c.score_type ) 
end

#
# Function that returns the scores written in the compactlog file 
# of one of the models given its index
#

function scores( c :: CompactLog, index :: Int64 )
  x = Vector{Float64}(undef,c.nmodels) 
  for jndex in 1:c.nmodels
    ipos = icl(c,index,jndex)
    if ipos == 0
      x[jndex] = 1.
      continue
    end
    x[jndex] = c.score[ipos]
  end
  return x
end

function scores( c :: CompactLog, model :: String )
  index = findfirst(isequal(model),c.name)
  return scores( c, index )
end

#
# Function that returns the index in vector score of the compactlog
# matrix, given the indexes of the pair to be studied
#

function icl( c :: CompactLog, index :: Int64, jndex :: Int64)

  #upper triangular, without diagonal
  #i,j ->   (i-1)*n+j - (i*(i-i))/2 + i

  if index == jndex
    return 0
  elseif index > jndex
    j = index
    i = jndex
  elseif index < jndex
    i = index
    j = jndex
  end
  return (i-1)*c.nmodels+j - (div((i*(i-1)),2)+i)

end


