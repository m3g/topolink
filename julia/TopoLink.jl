module TopoLink

  export davis, gscore

  # Data structures
  include("./src/structures.jl")

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

  # Functions to easy ploting data
  include("./src/simpleget.jl")

end









