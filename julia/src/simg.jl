#
# Function to produce a similarity vs. gscore plot 
#

function simg( c :: CompactLog, cscore :: Vector{Cscore}, model :: String )
  index = findfirst(isequal(model),c.name)
  return simg( c, cscore, index )
end
function simg( c :: CompactLog , cscore :: Vector{Cscore}, index :: Int64  )
  x = Vector{Float64}(undef,c.nmodels)
  y = Vector{Float64}(undef,c.nmodels)
  for j in 1:c.nmodels
    x[j] = cscore[j].gscore
    if j == index
      y[j] = 1.
    else
      ipos = icl(c,index,j) 
      y[j] = c.score[ipos] 
    end
  end
  return x, y
end
