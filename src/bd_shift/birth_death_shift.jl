## Wrapper function
export birth_death_shift
export plottree

@doc raw"""
    birth_death_shift(model, data)

Calculates average branch rates under the birth-death-shift model with a finite state space.

Example:

```julia
using Pesto

phy = readtree(Pesto.path("bears.tre")) 
ρ = 1.0  
data = make_SSEdata(phy, "", ρ; include_traits = false) 
λ = [0.1, 0.2] 
μ = [0.05, 0.15] 

η = 0.05 
model = SSEconstant(λ, μ, η)

res = birth_death_shift(model, data)
```
"""
function birth_death_shift(model, data; nshifts = true, shift_bayes_factor = true) 
    Ds, Fs = backwards_forwards_pass(model, data)

    Ss = ancestral_state_probabilities(data, Ds, Fs)
    rates = tree_rates(data, model, Fs, Ss)

    if nshifts
        #nshift = compute_nshifts(model, data, Ds, Ss; ape_order = false)
        nshift = state_shifts_simple(model, data, Ds, Fs)
        append!(nshift, 0.0)
        rates[!,"nshift"] = nshift
    end

    if shift_bayes_factor
        bf = posterior_prior_shift_odds(model,data)
        append!(bf, NaN)
        rates[!,"shift_bf"] = bf
        rates[!,"shift_bf_log"] = log10.(bf)
    end

    return(rates)
end
