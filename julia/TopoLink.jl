module TopoLink

  # Functions to read and manipulate the compactlog file
  include("./compactlog.jl")

  # Function to return a square score matrix from the log data
  include("./dmatrix.jl")

  # Functions to compute consensus scores
  include("./degree.jl")
  include("./davis.jl")
  include("./gscore.jl")

  include("./cscore.jl")

  # Functions to extract data in easy plotable ways
  include("./simg.jl")

end









