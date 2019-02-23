#
# Function that computes the correlation between two links
#
# https://en.wikipedia.org/wiki/Phi_coefficient
#

function correlation( link1 :: LinkData, link2 :: LinkData )

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
  n1x = n[1,1] + n[1,2]
  n2x = n[2,1] + n[2,2]
  nx1 = n[1,1] + n[2,1]
  nx2 = n[1,2] + n[2,2]

  correlation = (n[1,1]*n[2,2] - n[1,2]*n[2,1]) / sqrt( n1x*n2x*nx1*nx2 )

  return correlation

end

function foundlink( status :: Int64 ) 
  if ( status == 0 || status == 1 || status == 5 ) 
    return true
  else
    return false
  end
end
