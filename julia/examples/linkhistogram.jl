#
# Run with: julia davis_nxl.jl
#

push!(LOAD_PATH,"../")
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

using TopoLink, Plots

models = TopoLink.models( "./data/loglist.txt" )

link = linkdata( models, 5 )

result, nmodels = linkresults( link )

bar(result,nmodels,xrotation=60,title=linkname(link))
bar!(xlabel="Result",ylabel="Number of models")
savefig("./linkhistogram.png")
println("Created plot: ./linkhistogram.png")


