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
  # Function that returns the scores of one of the models given its index
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

  # Function that returns the index in vector score of the compactlog
  # matrix, given the indexes of the pair to be studied

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

  #
  # Function that computes the degree of a model
  #

  function degree( c :: CompactLog, index :: Int64 ; cutoff = 0.5 )

    degree = 0

    n = c.nmodels
    for jndex in 1:c.nmodels

      ipos = icl(c,index,jndex)
      if ipos == 0
        continue
      end

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
    index :: Int64
    gscore :: Float64
    degree :: Int64
  end

  function gscore( c :: CompactLog ; cutoff = 0.5 )
    gscore = Vector{Gscore}(undef,c.nmodels)
    for i in 1:c.nmodels
      degree = TopoLink.degree( c, i, cutoff = cutoff )
      name = c.name[i]
      gscore[i] = Gscore(name, i, degree / ( c.nmodels-1 ), degree )
    end
    #sort!(gscore, by = x -> x.gscore, rev=true)
    return gscore
  end

  #
  # Function to produce a similarity vs. gscore plot 
  #

  function simg( c :: CompactLog, gscore :: Vector{Gscore}, model :: String )
    index = findfirst(isequal(model),c.name)
    return simg( c, gscore, index )
  end

  function simg( c :: CompactLog , gscore :: Vector{Gscore}, index :: Int64  )

    x = Vector{Float64}(undef,c.nmodels)
    y = Vector{Float64}(undef,c.nmodels)

    for j in 1:c.nmodels
      x[j] = gscore[j].gscore
      if j == index
        y[j] = 1.
      else
        ipos = icl(c,index,j) 
        y[j] = c.score[ipos] 
      end
    end

    return x, y
  end

end









