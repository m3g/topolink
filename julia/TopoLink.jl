module TopoLink

  export names, davis, gscore, degree, nconsist, index
  export linkensemble
  export status, result
  export linkdata
  export linkresults

  # Data structures
  include("./src/structures.jl")

  # Functions to get status index or name from one another
  include("./src/status.jl")
  include("./src/result.jl")

  # Functions to read and manipulate the compactlog file
  include("./src/compactlog.jl")

  # Function to return a square score matrix from the log data
  include("./src/dmatrix.jl")

  # Functions to compute consensus scores
  include("./src/degree.jl")
  include("./src/davis.jl")
  include("./src/gscore.jl")

  # Functions to read topolink log files
  include("./src/readlog.jl")

  # Function to setup all data for a model
  include("./src/model.jl")

  # Function to compute the linkensemble plot
  include("./src/linkensemble.jl")
  
  # Functions to get link data for all models
  include("./src/linkdata.jl")
  include("./src/linkresults.jl")

  # Correlation functions
  include("./src/correlation.jl")
  include("./src/cmatrix.jl")
  export correlation, cmatrix

  # Functions to easy ploting data
  include("./src/simpleget.jl")

  # Functions to write some data in specific formats
  include("./src/write.jl")
  export linkname

end









