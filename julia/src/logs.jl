#
# Functions to read topolink log files
#

struct Residue
  name :: String
  chain :: String
  number :: Int64
  atom :: String
  sa :: Bool # Solvent accessibility
end

struct Link
  resid1 :: Residue
  resid2 :: Residue
  euclidean :: Float64
  topological :: Float64
  observed :: Bool
  dmin :: Float64
  dmax :: Float64
end

struct TopoLinkLog
  pdb :: String
  nlinks :: Int64
  link :: Vector{Link}
end

function readlog( filename :: String )

  local nlinks :: Int64
  local pdb :: String
  local link :: Vector{Link}
  local ra 

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

      link[ilink] = Link( resid1, resid2, euclidean, topological, observed, dmin, dmax )

    end
  end
  close(file)

  return TopoLinkLog( pdb, nlinks, link )

end
