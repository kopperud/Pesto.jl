# [Simple analysis](@id simple)

Here is an example of an analysis of branch-specific rates under the birth-death-shift model.

## Tree file

First, we load the necessary modules and read in the tree file.

```@setup simple
using Pesto

ρ = 0.635

include("../../src/primates.jl")
```
```julia simple
using Pesto

phy = readtree(Pesto.path("primates.tre"))
ρ = 0.635
primates = SSEdata(phy, ρ)
```

## Analysis
A simple analysis can be done like so:
```@example simple
model, rates = pesto(primates)
nothing # hide
```
To see how this analysis is set up, see the next section ([Extended analysis](@ref extended)).

## Tree plots
If we want to plot the results immediately, we will need to use the `Makie` dependency. `Makie` has several backends, but we will use `CairoMakie` for now, as it is designed to plot high-quality 2-dimensional figures.

```@example simple
using Makie, CairoMakie

treeplot(primates, rates)
```
If you want to instead plot the number of shifts, and for example use a different color gradient, you can enter the following
```@example simple
cmap = Makie.cgrad([:black, :red])
f = treeplot(primates, rates, "nshift"; cmap = cmap)
```
The figure can be saved using the following code
```julia
Makie.save("path/to/figure.pdf", f)
```
See section [Plot with ggtree](@ref ggtree) for more instructions on how to make more customized and publication-quality figures.

## Save to file
A summary of the results of the analysis is stored in the `rates` object, which is a data frame. We can print the first five rows to get an idea of it:
```@example simple
rates[1:5,:]
```
Each row in the data frame represents a branch in the tree, indexed by the edge index `edge`. 
Alternatively, you can index the rows using the `node` index. 
It is possible to save the data-frame as is to a file, using `using CSV; CSV.write("rates.csv", rates)`. We can also represent the output as an extended newick string:
```@example simple
newick(primates, rates)
```
If we want to save the newick string to a file, we can use the `writenewick` function
```julia
writenewick("primates_analysis.tre", primates, rates)
```
This tree file can be loaded in other programs such as `R` and can be plotted using standard packages like `ape` and `ggtree` (see section [Plot with ggtree](@ref ggtree)).