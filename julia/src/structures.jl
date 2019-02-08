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
  result :: String
end

struct TopoLinkLog
  pdb :: String
  model :: String
  nlinks :: Int64
  link :: Vector{Link}
  nconsist :: Int64
  nnotcons :: Int64
  nmissing :: Int64
end

import Base: ==

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

function ==( x :: Link, y :: Link ) 
  if ( x.resid1 == y.resid1 && x.resid2 == y.resid2 ) ||
     ( x.resid1 == y.resid2 && x.resid2 == y.resid1 ) ||
     ( x.resid2 == y.resid1 && x.resid1 == y.resid2 )
    return true
  else
    return false
  end
end


 












