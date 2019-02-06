#
# Function to return a square distance matrix given the compaclog score data
#

function inverse_minus_one( score )
  if score == 0 
    error(" In function inverse_minus_one (1/score-1) the score cannot be zero. ")
  end
  return 1. / score - 1.
end

function dmatrix( c :: CompactLog ; f = inverse_minus_one )
  dmatrix = Matrix{Float64}(undef,c.nmodels,c.nmodels) 
  k = 0
  for i in 1:c.nmodels-1
    for j in i+1:c.nmodels
      k = k + 1
      dmatrix[i,j] = f(c.score[k])
      dmatrix[j,i] = dmatrix[i,j]
    end
  end
  for i in 1:c.nmodels
    dmatrix[i,i] = f(i,i)
  end
  return dmatrix
end
