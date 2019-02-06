#
# Function that computes the degree of a model
#

function degree( c :: CompactLog, index :: Int64 ; cutoff = 0.5 )
  degree = 0
  for jndex in 1:c.nmodels
    ipos = icl(c,index,jndex)
    if ipos == 0
      continue
    end
    if c.score[ipos] > cutoff
      degree = degree + 1
    end
  end
  return degree
end

function degree( c :: CompactLog, model :: String ; cutoff = 0.5 )
  index = findfirst(isequal(model),c.name)
  return degree( c, index, cutoff=cutoff )
end
