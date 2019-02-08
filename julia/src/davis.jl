#
# Function that computes Davis consensus score of a model
#

function davis( c :: CompactLog, index :: Int64 ; f = inverse_minus_one )
  davis = 0.
  for jndex in 1:c.nmodels
    ipos = icl(c,index,jndex)
    if ipos == 0
      continue
    end
    davis = davis + f(c.score[ipos])
  end
  davis = davis / (c.nmodels-1)
  return davis
end

function davis( c :: CompactLog, model :: String ; f = inverse_minus_one )
  index = findfirst(isequal(model),c.name)
  return davis( c, index, f = f )
end
