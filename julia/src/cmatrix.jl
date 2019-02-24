#
# Function that computes the correlation matrix between links
#
# https://en.wikipedia.org/wiki/Phi_coefficient
#

function cmatrix( models :: Vector{Model} ; get = "phi" )

  nmodels = length(models)
  nlinks = models[1].nlinks

  if get == "phi" 
    cmatrix = Matrix{Float64}(undef,nlinks,nlinks)
  else
    cmatrix = Matrix{Int64}(undef,nlinks,nlinks)
  end

  for i in 1:nlinks-1
    linki = linkdata( models, i )
    for j in i+1:nlinks
      linkj = linkdata( models, j )
      cmatrix[i,j] = correlation( linki, linkj, get = get )
      cmatrix[j,i] = cmatrix[i,j]
    end
  end

  if get == "phi" 
    for i in 1:nlinks
      cmatrix[i,i] = 1.
    end
  elseif get == "n11"
    for i in 1:nlinks
      linki = linkdata( models, i )
      result, n = linkresults( linki )
      cmatrix[i,i] = n[1] + n[2] + n[6]
    end
  else
    for i in 1:nlinks
      cmatrix[i,i] = 0
    end
  end
  return cmatrix  

end
