
# TopoLink.jl :  Tools to analyze ensembles of models with TopoLink

This is a dynamic set of tools to analyze ensembles of structural models studied with TopoLink

Requisites: the Julia language and the packages Plots and ProgressMeter.

How to use it:

Load the TopoLink julia package using:
```
push!(LOAD_PATH,"/path_to/topolink/julia")
using TopoLink
using Plots
```

After obtaining the compactlog file as explained in the main TopoLink
site, load it using:

```
compactlog = TopoLink.compactlog("./compactlog-TMscore.dat")
```

The TopoLink load data can be read into a vector using: 

```
models = TopoLink.models( "./loglist.txt", compactlog=compactlog )
```

where `loglist.txt` is a file containing the list of TopoLink log files generated.

# Model data

The `models` vector contains the data for all models, meaning, for example:

`model[1].name` : model name

`model[1].davis` : model Davis consensus score

`model[1].nconsist` : number of consistent links of this model

`model[1].degree` : model clustering degree

The vector can be sorted by any of these properties by using the standard Julia sort function:

`sort!( models, by = model -> model.degree, rev = true)`

If you are not familiar to Julia, this means that the vector `models`
will be sorted by a property consisting in a function which, given each
element of the vector, returns the `degree` property of that element.
`rev = true` is used here so that the models are sorted from greater to
smaller clustering degrees. After this sorting, `model[1].name` will be
the name of the model of greater degree.

Once the compactlog and model data are loaded, you can easily plot many
properties of the models as function of other properties. For instance,
to plot the number of links consistent with each model as a function of
their davis consensus score, just do:

```
using Plots
scatter( davis(models), nconsist(models) )
scatter!( xlabel="Davis score", ylabel="Number of consistent XLs")
```

This will produce the following plot:

<p align="center">
<img src="https://github.com/mcubeg/topolink/blob/master/julia/examples/davis_nxl.png?raw=true">
</p>

## Evaluating the number of models necessary to satisfy the crosslinks

We also provide the function necessary to compute the number of models,
as classified by any of the model features, necessary to cumulatively
satisfy the crosslinks. This is the `linkensemble` function:

```
nsatisfied = linkensemble(models, by = model -> model.degree)
```

and these data can be plotted as a function of the index of the model. The
data is sorted by the `linkensemble` function by the chosen model
property. To plot the result, do  

```
plot( index(models), nsatisfied, , linewidth=2 ) 
plot!(xlabel="Models ordered by degree", ylabel="Cumulative XL satisfied")
plot!(xlim=[0,50])
```

producing the linkensemble plot:

<p align="center">
<img src="https://github.com/mcubeg/topolink/blob/master/julia/examples/ensemble.png?raw=true">
</p>

To view the results of the modeling for a single link, you can use the
`linkdata` structure and the `linkresult` function. For example:

```
link = linkdata( models, 5 )
```
will get all the data concerning link number 5 into the `link`
structure. The results obtained can be summarized with the `linkresult`
function, as in
```
result, nmodels = linkresults( link )
```
This will return, the name of the result ("OK: FOUND", etc) and the
number of models for which that result was observed. 
An histogram with this data can be easily constructed with:
```
bar(result,nmodels,xrotation=60,title=linkname(link))
bar!(xlabel="Result",ylabel="Number of models")
```
resulting in the figure below:
<p align="center">
<img src="https://github.com/mcubeg/topolink/blob/master/julia/examples/linkhistogram.png?raw=true">
</p>

## Correlation between crosslinks









