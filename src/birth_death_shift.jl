## Wrapper function
export birth_death_shift
export plottree
export bds

@doc raw"""
    birth_death_shift(model, data[; verbose = false])

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
function birth_death_shift(model, data; verbose = false)
    Ds, Fs = backwards_forwards_pass(model, data; verbose = verbose)
    Ps = ancestral_state_probabilities(data, model, Ds, Fs)

    res = calculate_tree_rates(data, model, Ds, Fs, Ps)

    out = Dict()

    out["lambda"] = res["average_node_rates"]["λ"]
    out["mu"] = res["average_node_rates"]["μ"]

    return(out)
end

"""
    bds(model, data[; verbose = false])

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

res = bds(model, data)
```
"""
function bds(model, data; verbose = false)
    Ds, Fs = backwards_forwards_pass(model, data; verbose = verbose)
    Ps = ancestral_state_probabilities(data, model, Ds, Fs)

    res = calculate_tree_rates(data, model, Ds, Fs, Ps)


    lambda = res["average_node_rates"]["λ"]
    mu = res["average_node_rates"]["μ"]

    phy = Dict("edge" => data.edges,
      "tip.label" => data.tiplab,
      "Nnode" => length(data.tiplab)-1,
     "edge.length" => data.branch_lengths)

    out = SSEresult(phy, lambda, mu)

    return(out)
end

"""
    plottree(x) 

Example:

```julia
res = bds(model, data)
plottree(res)
```
"""
function plottree(x)
    phy = x.phy
    lambda = x.lambda
    mu = x.mu

    RCall.@rput lambda
    RCall.@rput mu
    RCall.@rput phy

    RCall.R"""
    class(phy) <- "phylo"
    th <- max(ape::node.depth.edgelength(phy))
    lambda_average <- mean(lambda)

    df1 <- tibble::tibble("node" = 1:max(phy$edge),
            "Speciation rate" = lambda_average)
    x <- tidytree::as_tibble(phy)

    phydf <- merge(x, df1, by = "node")
    td_phy <- tidytree::as.treedata(phydf)

    p1a <- ggtree::ggtree(td_phy, ggplot2::aes(color = `Speciation rate`)) +
        ggtree::geom_tiplab(size = 8) +
        ggplot2::theme(legend.position = c(0.2, 0.8)) +
        ggplot2::xlim(c(0.0, th + 10)) 
    plot(p1a)
    """
    return 0;
end




