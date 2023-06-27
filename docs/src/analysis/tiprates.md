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

Let's use the LogNormal quantiles again.

```@example tips
λml, μml = estimate_constant_bdp(primates)

H = 0.587
n = 6

using Distributions
dλ = LogNormal(log(λml), H)
dμ = LogNormal(log(µml), H)

λquantiles = make_quantiles(dλ, n)
µquantiles = make_quantiles(dμ, n)
λ, μ = allpairwise(λquantiles, µquantiles)
η = optimize_eta(λ, µ, primates)
model = SSEconstant(λ, μ, η)
nothing # hide
```



## Ancestral state probabilities

We calculate the expression for the ancestral state probabilities (`S(t)`).
```@example tips
Ds, Fs = backwards_forwards_pass(model, primates)
Ss = ancestral_state_probabilities(primates, Ds, Fs)
nothing # hide
```

The `Ss` is a dictionary indexed over the branches. Each item corresponds to the function that when evaluated at time `t` yields the ancestral state probabilities for being in each state. For example, at branch `i=5`, at time `t=3.1`, the ancestral state probabilities are
```@example tips
Ss[5](3.1)
```

## Tip rates

If we loop over the branch incides, and find out which branches are terminal branches, we can store the mean rate value at time `t=0` for each species:
```@example tips
branch_indices = 1:size(primates.edges)[1]
ntips = length(primates.tiplab)

tip_rates = Dict{String,Float64}()
for i in branch_indices
    parent, child = primates.edges[i,:]

    if child <= ntips
        tiplab = primates.tiplab[child]

        rate = model.λ
        t = Ds[i].t[1]
        S = Ss[i](t)
        tip_rates[tiplab] = rate' * S
    end
end
tip_rates
```

## Distribution

If we plot the tip rates as a histogram, we can see that the primates tips are bimodally distributed. The high-rate species are the Old World Monkeys, and the low-rate species is everything else in the tree.

```@example tips
histogram(collect(values(tip_rates)), bins = 30, 
          grid = false, label = "",
          xlabel = "Tip rate (speciation rate, λ)",
          ylabel = "Frequency", size = (500, 300))
```