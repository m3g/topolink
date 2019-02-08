import Base: ==

#
# Amino acid residue 
#

struct Residue
  name :: String
  chain :: String
  number :: Int64
  sa :: Bool # Solvent accessibility
end

function ==( x :: Residue, y :: Residue )
  if x.name == y.name &&
     x.chain == y.chain &&
     x.number == y.number &&
     x.atom == y.atom 
    return true
  else
    return false
  end
end

#
# Link Atom
#

struct LinkAtom
  name :: String 
  residue :: Residue
end 

function Base.show( io :: IO, atom :: LinkAtom )
  print( atom.residue.name," ",atom.residue.chain," ",atom.residue.number," ",atom.name )
end

#
# A specific link
#

struct Link
  atom1 :: LinkAtom
  atom2 :: LinkAtom
  euclidean :: Float64
  topological :: Float64
  observed :: Bool
  dmin :: Float64
  dmax :: Float64
  result :: String
end

function ==( x :: Link, y :: Link ) 
  if ( x.atom1 == y.atom1 && x.atom2 == y.atom2 ) ||
     ( x.atom1 == y.atom2 && x.atom2 == y.atom1 ) ||
     ( x.atom2 == y.atom1 && x.atom1 == y.atom2 )
    return true
  else
    return false
  end
end

#
# A topolink log file
#

struct TopoLinkLog

  pdb :: String
  name :: String
  nlinks :: Int64
  link :: Vector{Link}
  nconsist :: Int64
  nnotcons :: Int64
  nmissing :: Int64

end

#
# All data for a model
#

struct Model

  name :: String
  pdb :: String
  nlinks :: Int64

  log :: String
  link :: Vector{Link}
  nconsist :: Int64
  nnotcons :: Int64
  nmissing :: Int64

  gscore :: Float64
  degree :: Int64
  davis :: Float64

end



 












