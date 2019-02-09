
include("./structures.jl")

function linkensemble(models)
  
  nlinks = models[1].nlinks

  nsatisfied = [ 0 for i in 1:length(models) ]

  satisfied = [ false for i in 1:nlinks ]

  for imodel in 1:length(models)
    for ilink in 1:nlinks
      if models[imodel].link[ilink].result == "OK: FOUND"
        satisfied[ilink] = true
      end
    end
    nsatisfied[imodel] = count(satisfied)
  end

  return nsatisfied

end
