module TopoLink

  using ProgressMeter

  struct CompactLog
  
    file :: String
    nmodels :: Int64
    name :: Vector{String}
    score :: Vector{Float64}
  
    alignment_file :: String
    pdb_list :: String
    score_type :: String
  
  end
  
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

  #
  # Function to return a square distance matrix given the compaclog score data
  #
  
  function simple_distance( i :: Int64, j :: Int64 ; score = 0. )
    if i == j 
      return 0.
    else
      if score == 0 
        error(" Pair score is zero for pair ", i," ",j)
      end
      return 1. / score - 1.
    end
  end

  function dmatrix( c :: CompactLog ; f = simple_distance )
    dmatrix = Matrix{Float64}(undef,c.nmodels,c.nmodels) 
    k = 0
    for i in 1:c.nmodels-1
      for j in i+1:c.nmodels
        k = k + 1
        dmatrix[i,j] = f(i,j,score=c.score[k])
        dmatrix[j,i] = dmatrix[i,j]
      end
    end
    for i in 1:c.nmodels
      dmatrix[i,i] = f(i,i)
    end
    return dmatrix
  end

  #
  # Function that computes the degree of a model
  #

  function degree( c :: CompactLog, index :: Int64 ; cutoff = 0.5 )

    degree = 0

    #upper triangular, without diagonal
    #i,j ->   (i-1)*n+j - (i*(i-i))/2 + i

    n = c.nmodels
    for jndex in 1:c.nmodels

      if index == jndex
        continue
      elseif index > jndex
        j = index
        i = jndex
      elseif index < jndex
        i = index
        j = jndex
      end

      ipos = (i-1)*n+j - (div((i*(i-1)),2)+i)
      if c.score[ipos] > cutoff
        degree = degree + 1
      end

    end
    return degree

  end

  function degree( c :: CompactLog, model :: String ; cutoff = 0.5 )
    index = findfirst(isequal(model),c.name)
    return degree( c, index, cutoff=cutoff )
  end

  #
  # Function that returns all degrees and G-scores of all models
  #

  struct Gscore
    name :: String
    gscore :: Float64
    degree :: Int64
  end

  function gscore( c :: CompactLog ; cutoff = 0.5 )
    gscore = Vector{Gscore}(undef,c.nmodels)
    for i in 1:c.nmodels
      degree = TopoLink.degree( c, i, cutoff = cutoff )
      name = c.name[i]
      gscore[i] = Gscore(name, degree / ( c.nmodels-1 ), degree )
    end
    sort!(gscore, by = x -> x.gscore, rev=true)
    return gscore
  end

end









