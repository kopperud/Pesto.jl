# [Number of rate shifts](@id shifts)

In this vignette we go a little more in-depth and explain how the number of rate shifts ($\hat{N}$) is estimated. 

### Tree file

First, we load the necessary modules and read in the tree file.

```@setup shift
using Pesto
using CairoMakie

ρ = 0.635

include("../../src/primates.jl")
```
```julia shift
using Pesto
using CairoMakie

phy = readtree(Pesto.path("primates.tre"))
ρ = 0.635
primates = SSEdata(phy, ρ)
```

### Model setup

In this vignette, we pick the rate values by hand, and we don't use so many, in order to illustrate how the calculations work.

```@example shift
tree_length = sum(primates.branch_lengths)
λ = [0.1, 0.2, 0.3, 0.4, 0.20]
μ = [0.05, 0.15, 0.05, 0.15, 0.25]
η = 1 / tree_length
model = SSEconstant(λ, μ, η)
nothing # hide
```

## Number of total rate shifts

If we want to see the total number of rate shifts per branch, we can use the function `birth_death_shift` for the standard inference:

```@example shift
rates = birth_death_shift(model, primates)
rates[1:5,:]
```

We can use `Makie` to plot the number of accumulated diversification rate shifts on the branches of the phylogeny
```@example shift
using Makie, CairoMakie

cmap = Makie.cgrad([:black, :red])
treeplot(primates, rates, "nshift"; cmap = cmap)
```

## Number of rate shifts

The total number of rate shifts, as shown above, is a reduction or simplification of the number of rate shifts that can be inferred by `Pesto`.
The number of rate shifts from state `j` to state `i` accumulated over the branch length (from old to young) is described by the following differential equation
```math
\frac{d\hat{N}_{M,ij}}{dt} = S_{M,j}(t) \frac{-\eta}{K-1} \frac{D_{M,i}(t)}{D_{M,j}(t)} \text{ if } j \neq i
```
with initial condition $\hat{N}_{ij}(t_0) = 0$. In Pesto, we would compute this using
```@example shift
nshift = state_shifts(model, primates)
nothing; # hide
```
The object returned `nshift` is a three-dimensional array. The first dimension corresponds to the branch index (what was `M`). The second dimension represents the arrival state (`i`), and the third dimension represents the departure state (`j`). The order of the first index (`M`) is the same order as the edge indices in a standard `ape`-style tree in R.

If we sum over second and third dimension, we get the number of rate shifts per branch:
```@example shift
sum(nshift, dims = 2:3)[:,1,1]
``` 

If instead we sum over the first dimension, we get a breakdown over which rate transitions were more frequent:
```@example shift
Nmatrix = sum(nshift, dims = 1)[1,:,:]
``` 
In this case, the most frequent rate shift was from state `2` to state `4`, with $\hat{N} = 0.95$ number of rate shifts. Going from state `2` to state `4` under this model means an increase of $0.4-0.2=0.2$ in speciation rate units. This can for example be visualized using a histogram:
```@example shift
mids, bins = makebins(Nmatrix, model, -0.35, 0.35; nbins = 7)

f = Makie.Figure(size = (500, 300))
ax = Makie.Axis(f[1,1], 
        ylabel = "Number of rate shifts",
        xlabel = "Change in speciation rate (λi - λj)",
        xticks = round.(mids; digits = 2),
        xgridvisible = false,
        ygridvisible = false)
Makie.barplot!(ax, mids, bins[:,1], color = :maroon)
f
```
Most of the rate shift events represent a shift from a smaller to a larger speciation rate (i.e. $\lambda_i - \lambda_j > 0$), however some rate shifts are in the other direction ($\lambda_i - \lambda_j < 0$). There are also a few rate shift events where the speciation rate does not change ($\lambda_i - \lambda_j = 0$). In these events, it is the extinction rate that changes, and not the speciation rate. If we are interested in the question, "what is the overall magnitude of shifts?", we can calculate the mean shift magnitude (weighted by their frequencies):
```math
\frac{1}{\sum_{i,j}\hat{N}_{ij}} \sum_{i,j} \hat{N}_{ij} (\lambda_i - \lambda_j).
```
In Julia we can calculate it like so
```@example shift
shifts_weighted = Nmatrix .* (λ .- λ')
mean_magnitude = sum(shifts_weighted) / sum(Nmatrix)
```
meaning that the overall shift magnitude for the primates tree under this model was an increase of 0.098 speciation rate units.
