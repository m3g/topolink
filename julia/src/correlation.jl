#
# Function that computes the correlation between two links
#
# https://en.wikipedia.org/wiki/Phi_coefficient
#

function correlation( link1 :: LinkData, link2 :: LinkData ; get :: String = "phi" )

  nmodels = length(link1.status)
  n = zeros(2,2)

  for i in 1:nmodels
    s1 = link1.status[i]
    s2 = link2.status[i]
    if foundlink(s1) 
      if foundlink(s2)
        n[1,1] = n[1,1] + 1
      else
        n[1,2] = n[1,2] + 1
      end
    else
      if foundlink(s2)
        n[2,1] = n[2,1] + 1
      else
        n[2,2] = n[2,2] + 1
      end
    end
  end

  if get == "phi" 
    n1x = n[1,1] + n[1,2] + 1
    n2x = n[2,1] + n[2,2] + 1
    nx1 = n[1,1] + n[2,1] + 1
    nx2 = n[1,2] + n[2,2] + 1
    correlation = (n[1,1]*n[2,2] - n[1,2]*n[2,1]) / sqrt( n1x*n2x*nx1*nx2 )
    return correlation
  elseif get == "n11"
    return n[1,1]
  elseif get == "n00"
    return n[2,2]
  elseif get == "n10"
    return n[1,2]
  elseif get == "n01"
    return n[2,1]
  else
    error(" Correlation type must be: phi, n11, n00, n01, or n10. ")
  end
 
end

function foundlink( status :: Int64 ) 
  if ( status == 0 || status == 1 || status == 5 ) 
    return true
  else
    return false
  end
end
