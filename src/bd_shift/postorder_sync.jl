export postorder_sync

## synchronous log likelihood calculation
@doc raw"""
    postorder_async(data)

performs the postorder iteration synchronously or in series (left and right subtrees are run on the same thread)

Example:
```julia
using Pesto
phy = readtree(Pesto.path("primates.tre"))
sampling_probability = 0.635
data = make_SSEdata(phy, sampling_probability)

λ = [0.3, 0.15]
μ = [0.1, 0.2]
η = 0.01

model = BDSconstant(λ, μ, η)
E = extinction_probability(model, data)

D, sf = postorder_sync(model, data, E)
```
The `D` is the partial likelihood at the root node. `D` is re-scaled by a factor such that it sums to one. The rescaling factor is `sf`, which is given on a log scale

```julia
julia> D
2-element Vector{Float64}:
 0.013941190827498735
 0.9860588091725013

julia> sf
-705.9668193580866
```
"""
function postorder_sync(model::Model, data::SSEdata, E)
    ## Pre-compute descendants in hashtable
    descendants = make_descendants(data)

    K = number_of_states(model)

    ## Postorder traversal: computing the branch probabilities through time
    nrows = size(data.edges, 1)
    Ntip = length(data.tiplab)
    
    ## Storing the solution at the end of the branch
    elt = eltype(model)

    ## Storing the scaling factors
    pD = (model.λ, model.μ, model.η, K, E)
    u0 = ones(elt, K)
    tspan = (0.0, 1.0)
    ode = backward_prob(model)
    prob = OrdinaryDiffEq.ODEProblem{true}(ode, u0, tspan, pD)

    root_index = Ntip+1
    left_edge, right_edge = descendants[root_index]

    D_left, sf_left = subtree_sync(left_edge, prob, model, data, descendants, Ntip, K, elt)
    D_right, sf_right =  subtree_sync(right_edge, prob, model, data, descendants, Ntip, K, elt)

    D = D_left .* D_right
    sf = sf_left + sf_right

    c = sum(D)
    D = D ./ c
    sf = sf_left + sf_right

    if c > 0.0
        sf += log(c)
    else
        sf -= Inf 
    end

    return(D, sf)
end

function subtree_sync(edge_index, prob, model, data, descendants, Ntip, K, elt)
    anc, dec = data.edges[edge_index,:]
    alg = OrdinaryDiffEq.Tsit5()  
    node_age = data.node_depth[dec]
    parent_node_age = data.node_depth[anc]
    tspan = (node_age, parent_node_age)  

    local sf_left, sf_right
    if dec <= Ntip
        u0 = ones(elt, K) .* data.sampling_probability
        sf_left = 0.0
        sf_right = 0.0
    else
        left_edge, right_edge = descendants[dec]
        
        D_left, sf_left = subtree_sync(left_edge, prob, model, data, descendants, Ntip, K, elt)
        D_right, sf_right = subtree_sync(right_edge, prob, model, data, descendants, Ntip, K, elt)

       
        λt = get_speciation_rates(model, node_age) 
        u0 = D_left .* D_right .* λt
    end

    sf = sf_left + sf_right
    if isinf(sf) | isnan(sf) | in(u0, NaN)
        D = ones(elt, K)
        sf -= Inf
    else
        prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
        sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg, 
                                   save_everystep = false, reltol = 1e-3)
        D = sol.u[end]
        c = sum(D)
        D = D ./ c

        if c > 0.0
            sf += log(c)
        else
            sf -= Inf
        end

        if !(sol.retcode == OrdinaryDiffEq.ReturnCode.Success)
            sf -= Inf
        end
    end

    return(D, sf)
end

