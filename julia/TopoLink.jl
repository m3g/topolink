module TopoLink

  # Functions to read and manipulate the compactlog file
  include("./src/compactlog.jl")

  # Function to return a square score matrix from the log data
  include("./src/dmatrix.jl")

  # Functions to compute consensus scores
  include("./src/degree.jl")
  include("./src/davis.jl")
  include("./src/gscore.jl")

  include("./src/cscore.jl")

  # Functions to extract data in easy plotable ways
  include("./src/simg.jl")

  # Functions to read topolink log files
  include("./src/logs.jl")

end









