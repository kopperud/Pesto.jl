## Wrapper function
export birth_death_shift
export plottree

@doc raw"""
    birth_death_shift(model, data)

Calculates average branch rates under the birth-death-shift model with a discrete rate categories.

Example:

```julia
using Pesto

phy = readtree(Pesto.path("bears.tre")) 
sampling_probability = 1.0  
bears = SSEdata(phy, sampling_probability)

λ = [0.1, 0.2] 
μ = [0.05, 0.12] 
η = 0.05 
model = BDSconstant(λ, μ, η)

rates = birth_death_shift(model, bears)
```
"""
function birth_death_shift(model::Model, data::SSEdata; nshifts = true, shift_bayes_factor = true) 
    Ds, Fs = backwards_forwards_pass(model, data)

    Ss = ancestral_state_probabilities(data, Ds, Fs)
    rates = tree_rates(data, model, Fs, Ss)

    if nshifts
        #nshift = compute_nshifts(model, data, Ds, Ss; ape_order = false)
        nshift = state_shifts_simple(model, data, Ds, Fs)
        append!(nshift, 0.0)
        rates[!,:nshift] = nshift
    end

    if shift_bayes_factor
        bf = posterior_prior_shift_odds(model,data)
        append!(bf, NaN)
        rates[!,"shift_bf"] = bf
        rates[!,"shift_bf_log"] = log10.(bf)
    end

    return(rates)
end


function birth_death_shift(model::Model, tree::Root; nshifts = true, shift_bayes_factor = true) 
    Ds, Fs = backwards_forwards_pass(model, tree)
    Ss = ancestral_state_probabilities(tree, Ds, Fs)
    rates = tree_rates(tree, model, Fs, Ss)

    if nshifts
        #nshift = compute_nshifts(model, data, Ds, Ss; ape_order = false)
        nshift = state_shifts_simple(model, tree, Ds, Fs)
        append!(nshift, 0.0)
        rates[!,:nshift] = nshift
    end

    if shift_bayes_factor
        bf = posterior_prior_shift_odds(model,tree)
        append!(bf, NaN)
        rates[!,"shift_bf"] = bf
        rates[!,"shift_bf_log"] = log10.(bf)
    end

    return(rates)
end

