# [Extended rate analysis](@id extended)

In this vignette, we will do the same as in the simple analysis, but we explain how the model is set up in more detail.

## Tree file

First, we load the necessary modules and read in the tree file. We assume that we know the total number of extant species, and using that we can calculate the sampling fraction.

```@setup extended
using Pesto

sampling_fraction = 0.635

include("../../src/primates.jl")
```
```julia extended
using Pesto

phy = readtree(Pesto.path("primates.tre"))
sampling_fraction = 0.635
primates = SSEdata(phy, sampling_fraction)
```

## Model choice

Next, we set up the SSE model, including its dimensionality and hyperparameters. For this model, we will draw the speciation rate (λ) and extinction rate (µ) from LogNormal distributions. We pick the median of the LogNormal distributions such that they correspond to the maximum-likelihood estimates of the constant-rate birth-death model. This is also called the "empirical Bayes" approach, where we use the data twice. We pick the log-sd as `H = 0.587`, which corresponds to a LogNormal distribution whose 2.5%-97.5% quantile spans one order of magnitude. 

```@example extended
λml, μml = estimate_constant_bdp(primates)

H = 0.587
n = 6

using Distributions
dλ = LogNormal(log(λml), H)
dμ = LogNormal(log(µml), H)

λquantiles = make_quantiles(dλ, n)
µquantiles = make_quantiles(dμ, n)
λ, μ = allpairwise(λquantiles, µquantiles)
nothing # hide
```
The scatter plot of `λ` on the x-axis, and `µ` on the y-axis looks like the figure below (blue dots), with the quantiles of the LogNormal distributions on the margin.

![primatestree](../assets/quantiles.svg)

Next, we estimate the rate shift parameter η under the SSE model, conditional on λ and µ.
```@example extended
η = optimize_eta(λ, µ, primates)
```

The units of $\eta$ are number of rate shift events per lineage per time. The product of the tree length (the sum of all branch lengths) times $\eta$ will give us the number of expected rate shifts under the prior:
```@example extended
sum(primates.branch_lengths) * η
```

This allows us to set up the SSE model object:
```@example extended
model = SSEconstant(λ, μ, η)
nothing # hide
```

With the model and data objects we can for example calculate the log likelihood
```@example extended
logL_root(model, primates)
```

## Branch rates and shifts
Or we can compute both the postorder and preorder pass, and get the expected speciation and extinction rates per branch. The result is a data frame object, and we print the first five rows:
```@example extended
rates = birth_death_shift(model, primates)
rates[1:5,:]
```

## Tree plots

As before, we can use `Makie` to make some quick tree plots. Here we are plotting the average net-diversification rate per branch, with a two-color scheme going from black to green.
```@example extended
using Makie, CairoMakie

cmap = Makie.cgrad([:black, :green])
treeplot(primates, rates, "mean_netdiv"; cmap = cmap)
```
