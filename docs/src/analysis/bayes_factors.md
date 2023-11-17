# [Bayes factors](@id bayesfactor)

In this vignette we will assess the support for that a branch has at least one rate shift, versus the null hypothesis of no shifts. We will test this using Bayes factors, i.e. the relative support for one or more versus no shifts.

### Tree file

First, we load the necessary modules and read in the tree file.

```@setup bayes
using Pesto
using Plots

ρ = 0.635

include("../../src/primates.jl")
```
```julia bayes
using Pesto
using Plots

phy = readtree(Pesto.path("primates.tre"))
ρ = 0.635
primates = SSEdata(phy, ρ)
```

### Model setup

Let's use the standard `Pesto` analysis with empirical Bayes for the speciation and extinction rates, and maximum-likelihood shift rate.

```@example bayes
model, rates = pesto(primates; n = 6)
nothing # hide
```

## Plotting Bayes factors

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

We can for example plot the Bayes factor directly on the tree. Since the Bayes factor can vary considerably, we instead plot the log-transformed Bayes factors, which are more concentrated around 0. A log Bayes factor with a value of 0 means that the prior and posterior support for the shift hypotheses are equal. A value much larger than 0 means that there is support for one or more shifts. A value of less than 0 means that there is more support for 0 rate shifts.

```julia
R"""
library(ggplot2)
library(ggtree)
p1 <- ggtree(td, aes(color = log(shift_bf))) +
    scale_color_gradient2(low = "white", mid = "black", high = "red", midpoint = 0) +
    geom_tiplab(size=2) +
    labs(color = "log Bayes factor")
"""
```

Alternatively, we can assess which branches had a strong support for there being at least one shift. If we choose a significance level, for example that the Bayes factor has to be at least 10, we can compute for which branches that is the case.
R"""
significance_level <- 10

td@data$signif_shifts <- factor(td@data$shift_bf > significance_level)
p2 <- ggtree(td, aes(color = signif_shifts)) +
    scale_color_manual(values = c("black", "red"), labels = c("no support", "strong support")) +
    geom_tiplab(size=2) +
    labs(color = "Significant shift")
"""








