#
# Function that computes the correlation matrix between links
#
# https://en.wikipedia.org/wiki/Phi_coefficient
#

function cmatrix( models :: Vector{Model} )

  nmodels = length(models)
  nlinks = models[1].nlinks

  cmatrix = Matrix{Float64}(undef,nlinks,nlinks)

  for i in 1:nlinks-1
    cmatrix[i,i] = 1.
    linki = linkdata( models, i )
    for j in i+1:nlinks
      linkj = linkdata( models, j )
      cmatrix[i,j] = correlation( linki, linkj )
      cmatrix[j,i] = cmatrix[i,j]
    end
  end
  cmatrix[nlinks,nlinks] = 1.

  return cmatrix  

end
