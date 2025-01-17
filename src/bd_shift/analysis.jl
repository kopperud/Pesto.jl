## wrapper for analysis
export pesto
export pesto_twostep
export pesto_fossil

export empirical_bayes

function empirical_bayes(data::SSEdata; n = 6, sd = 0.587)
    λml, μml = estimate_constant_bdp(data)

    dλ = Distributions.LogNormal(log(λml), sd)
    dμ = Distributions.LogNormal(log(μml), sd)

    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)
    η = optimize_eta(λ, µ, data)
    model = BDSconstant(λ, μ, η)
    return(model)
end

@doc raw"""
pesto_twostep(data)

this function performs a two-step analysis, where the estimation of the
parameters is split in two steps:

### step one
* we set up two variables (``\hat{\lambda}`` and ``\hat{\mu}``) representing the mean birth and death rates, and we find an estimate for them by maximizing the likelihood of the parameters under the lineage-homogeneous birth-death process

### step two
* we set up a birth-death-shift model where the means of the base distributions for the speciation and extinction rates are ``\hat{\lambda}`` and ``\hat{\mu}``, and the shift rate (η) is a parameter
* we find the estimate for η by maximizing the likelihood under the birth-death-shift process

### branch results
finally, we compute branch-specific results, including i) mean rates, ii) mean no. rate shifts, iii) Bayes factors, given the model found by the two-step procedure



Example:
```julia
using Pesto

phy = readtree(Pesto.path("primates.tre"))
sampling_probability = 0.635
data = SSEdata(phy, sampling_probability)

model, rates = pesto_twostep(data)
```
"""
function pesto_twostep(data::SSEdata; n = 6, sd = 0.587)
    model = empirical_bayes(data; n = n)

    rates = birth_death_shift(model, data);
    return(model, rates)
end



@doc raw"""
pesto(data)

performs maximum-likelihood estimation of the base distribution
parameter and the shift rate, and calculates branch-specific
posterior mean i) rates, ii) no. shifts and iii) Bayes factors

Example:
```julia
using Pesto

phy = readtree(Pesto.path("primates.tre"))
sampling_probability = 0.635
data = SSEdata(phy, sampling_probability)

model, rates = pesto(data)
```
"""
function pesto(data::SSEdata; n = 6, sd = 0.587)
    optres, model, i = fit_BhDh(data; n = n, sd = sd, n_attempts = 5)

    rates = birth_death_shift(model, data);
    return(model, rates)
end

@doc raw"""
pesto_fossil(tree)

similar to `pesto(...)`, however for a fossilized birth-death-shift (FBDS) process.

This function does the following
1. sets up a FBDS model where the speciation rate and fossilization rate is allowed to vary, but extinction rate is assumed to be constant
2. estimates five parameters using maximum likelihood: i) the mean speciation rate, ii) the mean fossilization rate, iii) the (constant) extinction rate), iv) the speciation shift rate, and v) the fossilization shift rate. 
3. compute branch-specific results given the model with the ML estimates


Example:
```julia
using Pesto

phy = readtree(Pesto.path("osteoglossomorpha.tre"))
sampling_probability = 0.234
tree = construct_tree(phy, sampling_probability)

## replace extant leaves that are at least 0.5 Ma older
## than the present with terminal fossil sampling events
## (with survival probability = 1)
assign_fossils!(tree, 0.5)

model, rates = pesto_fossil(tree)
```
"""
function pesto_fossil(tree::Root; n = 6, sd = 0.587)
    optres, model, i = optimize_hyperparameters2(tree; n = n, sd = sd, n_attempts = 5)

    #rates = birth_death_shift(model, data);
    rates = tree_rates(tree, model); 
    return(model, rates)
end
