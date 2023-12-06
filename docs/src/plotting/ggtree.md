# [Plot with ggtree](@id ggtree)

Here is an example of how the results can be plotted using `ggtree` in `R`.

## Tree file

First, we load the necessary modules, read in the tree file, and do a quick analysis.

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
model, rates = pesto(primates)
nothing # hide
```

## Tree plots in Makie
If we want to plot the results immediately, we will need to use the `Makie` dependency. `Makie` has several backends, but we will use `CairoMakie` for now, as it is designed to plot high-quality 2-dimensional figures.

```@example simple
using Makie, CairoMakie

treeplot(primates, rates)
```

## Tree plots in ggtree

The plotting functionality include in `Pesto` is intended as an interactive tool, and it's functionality is rather rudimentary. If you wish to make a publication-quality plot, we recommend instead to use a more established, well-tested and feature-rich library especially for this purpose. A good choice is to use the `ggtree` library in `R`, although other libraries could be used as well.
In `Pesto`, we provide functionality to save the results as a Newick strinc with metadata for each node. If we want to save the newick string to a file, we can use the `writenewick` function
```julia
writenewick("primates_analysis.tre", primates, rates)
```
This tree file can be loaded in other programs such as `R` and can be opened using standard packages like `ape` and `treeio` or any other software that can handle extended Newick trees.

```R
library(treeio)
phy <- treeio::read.beast.newick("primates.analysis.tre")

library(ggplot2)
library(ggtree)
p1 <- ggtree(td, aes(color = mean_lambda)) +  
    geom_tiplab(size=2) +
    labs(color = "Mean speciation rate")
```

![primatestree](../assets/primates_lambda.svg)