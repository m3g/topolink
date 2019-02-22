# 
# Functions to read topolink log files
#

function readlog( filename :: String )

  local nlinks
  local pdb 
  local link
  local ra 
  local modelname
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
      modelname = name[1]
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
      residue = Residue( name, chain, number, sa )
      atom1 = LinkAtom( atom, residue )
     
      name = data[6]
      chain = data[7]
      number = parse(Int64,data[8])
      atom = data[9]
      if line[ra[2]:ra[2]] == "Y"
        sa = true
      else
        sa = false
      end
      residue = Residue( name, chain, number, sa )
      atom2 = LinkAtom( atom, residue )
   
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

      link[ilink] = Link( atom1, atom2, euclidean, topological, observed, dmin, dmax, result )

    end
  end
  close(file)

  return TopoLinkLog( pdb, modelname, nlinks, link, nconsist, nnotcons, nmissing )

end

function Base.show( io :: IO, m :: TopoLinkLog )
  print( m.name, " log file, with ", m.nlinks," links.")
end



