# [Tip rates](@id tiprates)

In this vignette we will calculate the diversification rates at the tips of the phylogeny.

### Tree file

First, we load the necessary modules and read in the tree file.

```@setup tips
using Pesto
using Plots

ρ = 0.635

include("../../src/primates.jl")
```
```julia tips
using Pesto
using Plots

phy = readtree(Pesto.path("primates.tre"))
ρ = 0.635
primates = SSEdata(phy, ρ)
```

### Model setup

Let's use the standard `Pesto.jl` model again, with LogNormal quantiles and maximum-likelihood shift rate ($\eta$).

```@example tips
model, rates = pesto(primates)
nothing # hide
```

## Tip rates

In order to calculate the value of the rates at the present, we need to evalute the following expression

```math
    \lambda_\text{tip} = \mathbf{S}(t=0)^\top \boldsymbol{\lambda},
```
where $\mathbf{S}(t=0)$ are the posterior state probabilities at time $t=0$, i.e. the present. We can compute the tip rates conveniently with the `tip_rates()` function, which gives a `DataFrame` as a result.

```@example tips
df = tip_rates(model, primates)
df[1:5,:]
```

## Distribution

If we plot the tip rates as a histogram, we can see that the primates tips are bimodally distributed. The high-rate species are the Old World Monkeys, and the low-rate species is everything else in the tree.

```@example tips
plots = []
for rate in [:lambda, :mu, :netdiv, :relext]
    p = histogram(df[!,rate], bins = 30, 
            grid = false, label = "",
            xlabel = string("Tip rate (", rate, ")"),
            ylabel = "Frequency", size = (500, 300))
    append!(plots, [p])
end
p3 = plot(plots...; layout = (2,2))
```

