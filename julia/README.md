
This will be a set of tools to analyze ensemble of models.

Currently, you can load all data for the models using:

```
include("./TopoLink.jl")`

compactlog = TopoLink.compactlog("./compactlog-TMscore.dat")

models = TopoLink.models( "./loglist.txt", compactlog=compactlog )
```

And plot, for example, the number of consistent links of each model as a function of their davis
consensus scores with:

```
using Plots
plot(TopoLink.davis(models),TopoLink.nconsist(models))

```
