#
# Run with: julia linkhistogram.jl
#

push!(LOAD_PATH,"../")
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

using TopoLink, Plots

# Read all logs into models vector
models = TopoLink.models( "./data/logs" )

link = linkdata( models, 5 )

result, nmodels = linkresults( link )

bar(result,nmodels,xrotation=60,label="")
bar!(xlabel="Result",ylabel="Number of models",title=linkname(link))
savefig("./linkhistogram.png")
println(" Created plot: ./linkhistogram.png")


