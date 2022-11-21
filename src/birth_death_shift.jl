## Wrapper function
export birth_death_shift

"""
    birth\\_death\\_shift(model, data[; verbose = false])

Calculates average branch rates under the birth-death-shift model with a finite state space.

Example:

using Diversification

treefile = "data/bears.tre"\\
phy = readtree(treefile) \\
ρ = 1.0  \\
data = make\\_SSEdata(phy, "", ρ; include_traits = false) \\
λ = [0.1, 0.2] \\
μ = [0.05, 0.15] \\

η = 0.05 \\
model = SSEconstant(λ, μ, η)

res = birth\\_death\\_shift(model, data)
"""
function birth_death_shift(model, data; verbose = false)
    Ds, Fs = backwards_forwards_pass(model, data; verbose = verbose)
    res = calculate_tree_rates(data, model, Ds, Fs; verbose = verbose)

    out = Dict()

    out["lambda"] = res["average_node_rates"]["λ"]
    out["mu"] = res["average_node_rates"]["μ"]

    return(out)
end
