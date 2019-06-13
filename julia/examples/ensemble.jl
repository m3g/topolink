#
# Run with: julia ensemble.jl
#

push!(LOAD_PATH,"../")
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

using TopoLink, Plots

compactlog = TopoLink.compactlog("./data/compactlog-TMscore.dat")

# Read all logs and alignment data into models vector
models = TopoLink.models( "./data/logs" , compactlog=compactlog )


# Compute the number of links satisfied by subsets of models ordered
# by degree. rev=true makes sense because the degree is greater for better models
nsatisfied = linkensemble(models, by = model -> model.degree, rev=true)

plot(index(models),nsatisfied,linewidth=2,label="") 
plot!(xlabel="Models ordered by degree", ylabel="Cumulative XL satisfied")
plot!(xlim=[0,50])
savefig("ensemble.png")

println(" Created plot: ensemble.png")


