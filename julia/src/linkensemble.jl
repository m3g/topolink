
include("./structures.jl")

function linkensemble( logs , cscore :: Cscore )
  
  nlinks = logs[1].nlinks
  nsatisfied = [ 0 for i in 1:length(logs) ]
  satisfied = [ false for i in 1:nlinks ]

  ilog = 0
  for log in logs
    ilog = ilog + 1
    for i in 1:log.nlinks
      if log.link[i] == log[1].link[i]

        if ! satisfied[log.

        if log.link[i].observed 
          if log.link[i].result == "OK: FOUND" 
            nsatisfied[ilog] = nsatisfied

          end
        end
      end
    end
  end

end
