#
# Function that returns the result string from the status integer
#

function result( status :: Int64 )

  if status == 0
    return "OK: FOUND"  
  end
  if status == 1
    return "BAD: SHORT"  
  end
  if status == 2
    return "BAD: LONG"  
  end
  if status == 3
    return "BAD: EUCL"  
  end
  if status == 4
    return "BAD: NOTFOUND"  
  end
  if status == 5
    return "BAD: MISSING"  
  end
  if status == 6
    return "OK: LONG"  
  end 
  if status == 7
    return "OK: EUCL"  
  end
  if status == 8
    return "OK: NOTFOUND"  
  end

end

