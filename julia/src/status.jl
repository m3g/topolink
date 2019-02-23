#
# Function that returns the status integer classifier given the result string
#

function status( result :: String )

  if result == "OK: FOUND" 
    return 0
  end
  if result == "BAD: SHORT" 
    return 1
  end
  if result == "BAD: LONG" 
    return 2
  end
  if result == "BAD: EUCL" 
    return 3
  end
  if result == "BAD: NOTFOUND" 
    return 4
  end
  if result == "BAD: MISSING" 
    return 5
  end
  if result == "OK: LONG" 
    return 6
  end 
  if result == "OK: EUCL" 
    return 7
  end
  if result == "OK: NOTFOUND" 
    return 8
  end

end

