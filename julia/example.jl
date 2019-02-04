
include("./TopoLink.jl")

file = "/home/leandro/Documents/fabio/salbIII/26rest/gscore/compactlog-TMscore.dat"

smatrix = TopoLink.compactlog(file) 

gscore = TopoLink.gscore(smatrix) 

for i in 1:10
  println(gscore[i].name," ",gscore[i].gscore)
end




