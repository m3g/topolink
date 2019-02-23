#
# Function that returns, for a link, the number of models for each type of result
#

function linkresults( linkdata :: LinkData )

  results = [ result(i) for i in 0:8 ]
  n = zeros(Int64,9) 
  for i in 1:9
    n[i] = count( status -> status == i-1, linkdata.status )
  end

  return results, n

end

