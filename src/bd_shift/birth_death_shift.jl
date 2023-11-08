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
        nshift = compute_nshifts(model, data, Ds, Ss; ape_order = false)
        append!(nshift, 0.0)
        rates[!,"nshift"] = nshift
    end

    if shift_bayes_factor
        bf = posterior_prior_shift_odds(model,data)
        append!(bf, NaN)
        rates[!,"shift_bf"] = bf
    end

    return(rates)
end

"""
    plottree(x) 

Example:

```julia
res = bds(model, data)
plottree(res)
```
"""
function plottree(data::SSEdata, rates::DataFrames.DataFrame)

    RCall.@rput rates
    p = data
    RCall.@rput p
    th = maximum(data.node_depth)
    RCall.@rput th

    RCall.R"""
    x <- tidytree::as_tibble(p)
    td <- tidytree::as.treedata(merge(x, rates, by = "node"))
    
    p1a <- ggtree::ggtree(td, ggplot2::aes(color = `mean_lambda`)) +
        ggtree::geom_tiplab(size = 8) +
        ggplot2::theme(legend.position = c(0.2, 0.8)) +
        ggplot2::xlim(c(0.0, th + 10)) 
    plot(p1a)
    """
    return 0;
end




