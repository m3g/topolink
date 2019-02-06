
include("../TopoLink.jl")

file = "/home/leandro/Drive/Work/topolink/compactlog-TMscore.dat"

smatrix = TopoLink.compactlog(file) 

consensus = TopoLink.cscores(smatrix) 

sort!(consensus, by = x -> x.gscore, rev=true)

for i in 1:10
  println(consensus[i].name," ",consensus[i].gscore)
end




