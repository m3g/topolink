#
# Run with: julia correlation.jl
#

push!(LOAD_PATH,"../")
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

using TopoLink, Plots

# Read all model data into models vector
models = TopoLink.models("./data/logs")

link1 = linkdata( models, 1 )
link2 = linkdata( models, 2 )

c12 = correlation( link1, link2 )

println(" Correlation between $(linkname(link1)) and $(linkname(link2)): ", c12)

C  = cmatrix( models )

heatmap(C,color=cgrad(:RdBu),clims=(-1,1))
heatmap!(xlabel="XL index",ylabel="XL index",title="XL correlation map")

savefig("./correlation.png")
println(" Created plot: ./correlation.png")

C  = cmatrix( models, get = "n11" )

heatmap(C,color=cgrad(:tempo))
heatmap!(xlabel="XL index",ylabel="XL index",title="Number of models satisfying both links")
savefig("./n11.png")
println(" Created plot: ./n11.png")


