
function davis( models :: Vector{Model} )
  return [ models[i].davis for i in 1:length(models) ]
end

function degree( models :: Vector{Model} )
  return [ models[i].degree for i in 1:length(models) ]
end

function gscore( models :: Vector{Model} )
  return [ models[i].gscore for i in 1:length(models) ]
end

function name( models :: Vector{Model} )
  return [ models[i].name for i in 1:length(models) ]
end

function nconsist( models :: Vector{Model} )
  return [ models[i].nconsist for i in 1:length(models) ]
end
