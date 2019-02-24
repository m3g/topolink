#
# Function that reads the data of a given link in all models and sets
# the linkdata vector with that data
#

function linkdata( models :: Vector{Model}, ilink :: Int64 )

  local result, status

  nmodels = length(models)
  nlinks = models[1].nlinks

  index = ilink
  atom1 = models[1].link[ilink].atom1
  atom2 = models[1].link[ilink].atom2
  dmin = models[1].link[ilink].dmin
  dmax = models[1].link[ilink].dmax
  observed = models[1].link[ilink].observed
  result = Vector{String}(undef,nmodels)
  status = Vector{Int64}(undef,nmodels)

  linkdata = LinkData(index,atom1,atom2,dmin,dmax,observed,result,status)

  for imodel in 1:nmodels
    linkdata.result[imodel] = models[imodel].link[ilink].result
    linkdata.status[imodel] = models[imodel].link[ilink].status
  end

  return linkdata

end

