#
# Run with: julia davis_nxl.jl
#

push!(LOAD_PATH,"../")
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

using TopoLink, Plots

models = TopoLink.models( "./data/loglist.txt" )

link1 = linkdata( models, 1 )
link2 = linkdata( models, 2 )

c12 = correlation( link1, link2 )

println(" Correlation between 1 and 2: ", c12)

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


