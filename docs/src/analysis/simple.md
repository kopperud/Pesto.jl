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
This tree file can be loaded in other programs such as `R` and can be plotted using standard packages like `ape` and `ggtree`.

## Tree plots
If we want to plot the results immediately, without saving it to a file, we can use the module `RCall` to access R packages directly. Julia objects can be exported to an R session using the macro `@rput`, (and retrieved from R with `@rget`). R code can be called by prefixing a string with `R`, e.g. `R"print()"`, or multiline `R"""..."""`. You can also enter the R session interactively through the Julia REPL by entering the character `$`. Here we plot the phylogeny using some R-packages that we load first.

```julia
using RCall

@rput primates
@rput rates

R"""
library(tibble)
library(tidytree)
x <- as_tibble(primates)
td <- as.treedata(merge(x, rates, by = "node"))
"""
```

We can plot the mean speciation rate

```julia
R"""
library(ggplot2)
library(ggtree)
p1 <- ggtree(td, aes(color = mean_lambda)) +  
    geom_tiplab(size=2) +
    labs(color = "Mean speciation rate")
"""
```
![primatestree](../assets/primates_lambda.svg)