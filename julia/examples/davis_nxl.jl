#
# Run with: julia davis_nxl.jl
#

push!(LOAD_PATH,"../")
ENV["GKSwstype"]="nul" # This supresses the need of a display while plotting

using TopoLink, Plots

compactlog = TopoLink.compactlog("./data/compactlog-TMscore.dat")
models = TopoLink.models( "./data/loglist.txt", compactlog=compactlog )

scatter( davis(models), nconsist(models) )
scatter!( xlabel="Davis score", ylabel="Number of consistent XLs")
savefig("davis_nxl.png")

println("Created plot: davis_nxl.png")


