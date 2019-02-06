#
# Structure and function to compute and contain consensus 
# measures scores
#

struct Cscore
  name :: String
  index :: Int64
  gscore :: Float64
  degree :: Int64
  davis :: Float64
end

function cscores( c :: CompactLog ; cutoff = 0.5, f = inverse_minus_one )
  cscore = Vector{Cscore}(undef,c.nmodels)
  for i in 1:c.nmodels
    degree = TopoLink.degree( c, i, cutoff = cutoff )
    davis = TopoLink.davis( c, i, f = f )
    gscore = degree / ( c.nmodels - 1 )
    name = c.name[i]
    cscore[i] = Cscore(name, i, gscore, degree, davis )
  end
  #sort!(cscore, by = x -> x.gscore, rev=true)
  return cscore
end
