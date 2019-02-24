
## TopoLink.jl :  Tools to analyze ensembles of models with TopoLink

This is a dynamic set of tools to analyze ensembles of structural models studied with TopoLink

* [Installation and data loading](#installation-and-data-loading)
* [Model data](#model-data) 
* [Evaluating the number of models necessary to satisfy the crosslinks](#evaluating-the-number-of-models-necessary-to-satisfy-the-crosslinks) 
* [Evaluating the ensemble properties of a specific link](#evaluating-the-ensemble-properties-of-a-specific-link)
* [Correlation between crosslinks](#correlation-between-crosslinks)  

### Installation and data loading

Requisites: the <a href=https://julialang.org/downloads/ target="_blank">Julia</a> language and the packages Plots and ProgressMeter.
Therefore, after installing the Julia interpreter, do:

```
using Pkg
Pkg.add("Plots","ProgressMeter")
```

When starting a job, load the TopoLink julia package using:
```
push!(LOAD_PATH,"/path_to/topolink/julia")
using TopoLink, Plots
```

The TopoLink log data can be read into a vector using: 
```
models = TopoLink.models( "./data/logs" )
```
where `./data/logs` is directory containing the TopoLink log files generated.

After obtaining the compactlog alignment file as explained in the main TopoLink
site, load it using:

```
compactlog = TopoLink.compactlog("./compactlog-TMscore.dat")
```

And, in this case, the TopoLink log data can be read including the alignment data, with

```
models = TopoLink.models( "./data/logs", compactlog=compactlog )
```

### Model data

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
their Davis consensus score, just do:

```
scatter( davis(models), nconsist(models) )
scatter!( xlabel="Davis score", ylabel="Number of consistent XLs")
```

This will produce the following plot:

<p align="center">
<img src="https://github.com/mcubeg/topolink/blob/master/julia/examples/davis_nxl.png?raw=true">
<br><br>
<a href="https://github.com/mcubeg/topolink/blob/master/julia/examples/davis_nxl.jl">
Complete Example
</a>
</p>

### Evaluating the number of models necessary to satisfy the crosslinks

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
plot( index(models), nsatisfied, linewidth=2 ) 
plot!(xlabel="Models ordered by degree", ylabel="Cumulative XL satisfied")
plot!(xlim=[0,50])
```

producing the linkensemble plot:

<p align="center">
<img src="https://github.com/mcubeg/topolink/blob/master/julia/examples/ensemble.png?raw=true">
<a href="https://github.com/mcubeg/topolink/blob/master/julia/examples/ensemble.jl">
<br><br>
Complete Example
</a>
</p>

### Evaluating the ensemble properties of a specific link

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
<a href="https://github.com/mcubeg/topolink/blob/master/julia/examples/linkhistogram.jl">
<br><br>
Complete Example
</a>
</p>

### Correlation between crosslinks

The `linkdata` structure can also be used to compute the correlation
between crosslinks. For example,

``` 
link1 = linkdata( models, 1 )
link2 = linkdata( models, 2 )
```
define two structures with the ensembles properties of these two links.
The correlation of these two links can then be computed with the
`correlation` function:
```
c12 = correlation( link1, link2 )
```
Of course, this feature can be easily extended to compute the complete
correlation matrix between the crosslinks. This is implemented in the
`cmatrix` function, to be used as:
```
C  = cmatrix( models )
```
This matrix can be plotted as a heatmap, using:
```
heatmap(C,color=cgrad(:RdBu),clims=(-1,1))
heatmap!(xlabel="XL index",ylabel="XL index",title="XL correlation map")
```
leading to the correlation matrix plot:
<p align="center">
<img src="https://github.com/mcubeg/topolink/blob/master/julia/examples/correlation.png?raw=true">
</p>

The `correlation` and `cmatrix` functions also accept a `get` argument
which by default is `"phi"` and returns the correlation coefficient.
`get` might be set, however, to `n11`, `n00`, `n10`, or `n01`, and the
functions will return the number of models satisfying both links,
neither of them, only the first one, or only the second one,
respectively.  

For example, a matrix of the number of models satisfying both links of
each possible pair simultaneously can be obtained with
```
C  = cmatrix( models, get = "n11" )
```
and this can be plotted with
```
heatmap(C,color=cgrad(:tempo))
heatmap!(xlabel="XL index",ylabel="XL index",title="Number of models satisfying both links")
```
leading to the following figure:
<p align="center">
<img src="https://github.com/mcubeg/topolink/blob/master/julia/examples/n11.png?raw=true">
<br><br>
<a href="https://github.com/mcubeg/topolink/blob/master/julia/examples/correlation.jl">
Complete Example
</a>
</p>












