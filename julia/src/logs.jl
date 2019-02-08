# 
# Functions to read topolink log files
#

include("./structures.jl")

function readlog( filename :: String )

  local nlinks
  local pdb 
  local link
  local ra 
  local model
  local nconsist, nnotcons, nmissing

  file = open(filename,"r")
  nlinks = 0
  for line in eachline(file)
    data = split(line) 
    if length(data) < 1 
      continue
    end
    if data[1] == "LINK:"
      nlinks = nlinks + 1
    end
  end
  link = Vector{Link}(undef,nlinks)

  seekstart(file)
  ilink = 0
  for line in eachline(file)
    data = split(line) 
    if length(data) < 1 ; continue ; end
    if data[1] == "RESIDUE1"
      ra = findfirst("RA",line)
    end
    if length(data) < 2 ; continue ; end
    if data[1] == "PDB" && data[2] == "input"
      pdb = data[4]
      name = basename(pdb)
      name = split(name,".")
      model = name[1]
    end
    if data[1] == "RESULT0:"
      nconsist = parse(Int64,data[2])
    end
    if data[1] == "RESULT2:"
      nnotcons = parse(Int64,data[2])
    end
    if data[1] == "RESULT3:"
      nmissing = parse(Int64,data[2])
    end
    if data[1] == "LINK:"
      ilink = ilink + 1

      name = data[2]
      chain = data[3]
      number = parse(Int64,data[4])
      atom = data[5]
      if line[ra[1]:ra[1]] == "Y"
        sa = true
      else
        sa = false
      end
      resid1 = Residue( name, chain, number, atom, sa )
     
      name = data[6]
      chain = data[7]
      number = parse(Int64,data[8])
      atom = data[9]
      if line[ra[2]:ra[2]] == "Y"
        sa = true
      else
        sa = false
      end
      resid2 = Residue( name, chain, number, atom, sa )
   
      euclidean = parse(Float64,data[10])
      if data[11][1:1] == ">"
        topological = 0.
      else
        topological = parse(Float64,data[11])
      end

      if data[12] == "YES"
        observed = true
      else
        observed = false
      end

      dmin = parse(Float64,data[13])
      dmax = parse(Float64,data[14])

      result = "$(data[15]) $(data[16])"

      link[ilink] = Link( resid1, resid2, euclidean, topological, observed, dmin, dmax, result )

    end
  end
  close(file)

  return TopoLinkLog( pdb, model, nlinks, link, nconsist, nnotcons, nmissing )

end

function Base.show( io :: IO, m :: TopoLinkLog )
  print( m.model, " log file, with ", m.nlinks," links.")
end

#
# Function to read a lot of log files
#

using ProgressMeter
function readlogs( loglistname :: String )

  loglist = open(loglistname,"r")

  nlogs = 0
  for filename in eachline(loglist)
    nlogs = nlogs + 1
  end
  seekstart(loglist)

  logs = Vector{TopoLinkLog}(undef,nlogs)

  seekend(loglist) ; fileSize = position(loglist) ; seekstart(loglist)
  p = Progress(fileSize,1," Reading list of log files: ")

  ilog = 0
  for filename in eachline(loglist)
    update!(p,position(loglist))
    ilog = ilog + 1
    logs[ilog] = readlog(filename)
  end

  close(loglist)

  return logs 

end

function Base.show( io :: IO, ::MIME"text/plain", logs :: Array{TopoLinkLog} )
  print(" Log list with ", length(logs)," log files.")
end



