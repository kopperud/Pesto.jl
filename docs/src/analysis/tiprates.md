# [Tip rates](@id tiprates)

In this vignette we will calculate the diversification rates at the tips of the phylogeny.

### Tree file

First, we load the necessary modules and read in the tree file. We assume that we know the total number of extant species, and using that we can calculate the sampling fraction.

```@setup tips
using Pesto
using CairoMakie

sampling_fraction = 0.635

include("../../src/primates.jl")
```
```julia tips
using Pesto
using CairoMakie

phy = readtree(Pesto.path("primates.tre"))
sampling_fraction = 0.635
primates = SSEdata(phy, sampling_fraction)
```

### Model setup

Let's use the standard `Pesto.jl` model again, with LogNormal quantiles and maximum-likelihood shift rate ($\eta$).

```@example tips
model, rates = pesto(primates)
nothing # hide
```

## Tip rates

In order to calculate the values of the rates at the present, we need to evalute the following expressions

```math
    \begin{aligned}
        \lambda_\text{tip} &= \mathbf{S}(t=0)^\top \boldsymbol{\lambda}\\
        \mu_\text{tip} &= \mathbf{S}(t=0)^\top \boldsymbol{\mu}\\
        r_\text{tip} &= \mathbf{S}(t=0)^\top (\boldsymbol{\lambda}-\boldsymbol{\mu})\\
        \epsilon_\text{tip} &= \mathbf{S}(t=0)^\top (\boldsymbol{\mu} \oslash \boldsymbol{\lambda}),
    \end{aligned}
```
where $\mathbf{S}(t=0)$ are the posterior state probabilities at time $t=0$, i.e. the present. We can compute the tip rates conveniently with the `tip_rates()` function, which gives a `DataFrame` as a result.

```@example tips
df = tip_rates(model, primates)
df[1:5,:]
```

## Distribution

If we plot the tip rates as a histogram, we can see that the primates tips are bimodally distributed. The high-rate species are the Old World Monkeys, and the low-rate species is everything else in the tree.

```@example tips
f = Figure(resolution = (300, 600))

axs = []
colors = [:blue, :red, :green, :orange]
for (i, rate) in enumerate([:lambda, :mu, :netdiv, :relext])
    ax = Axis(f[i,1], 
        xgridvisible = false,
        ygridvisible = false,
        xlabel = string("Tip rate (", rate, ")"),
        ylabel = "Frequency")

    hist!(ax, df[!,rate], bins=30, color = colors[i])
end
rowgap!(f.layout, 5.0)
f
```

