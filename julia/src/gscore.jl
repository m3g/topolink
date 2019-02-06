#
# Function that computes the degree of a model
#

function gscore( c :: CompactLog, index :: Int64 ; cutoff = 0.5 )
  gscore = 0.
  for jndex in 1:c.nmodels
    ipos = icl(c,index,jndex)
    if ipos == 0
      continue
    end
    if c.score[ipos] > cutoff
      gscore = gscore + 1.
    end
  end
  return gscore
end

function gscore( c :: CompactLog, model :: String ; cutoff = 0.5 )
  index = findfirst(isequal(model),c.name)
  return gscore( c, index, cutoff=cutoff )
end
