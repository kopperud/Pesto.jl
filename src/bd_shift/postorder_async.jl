export postorder_async

## async log likelihood calculation
@doc raw"""
    postorder_async(data)

performs the postorder iteration asynchronously (left and right subtrees are run on different threads)

Example:
```julia
using Pesto
phy = readtree(Pesto.path("primates.tre"))
sampling_probability = 0.635
data = make_SSEdata(phy, sampling_probability)

λ = [0.3, 0.15]
μ = [0.1, 0.2]
η = 0.01

model = SSEconstant(λ, μ, η)
E = extinction_probability(model, data)

D, sf = postorder_async(model, data, E)
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
function postorder_async(model::SSE, data::SSEdata, E)
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

    local sf_left, sf_right
    local D_left, D_right
    Threads.@sync begin
        Threads.@spawn D_left, sf_left = subtree(left_edge, prob, model, data, descendants, Ntip, K, elt)
        D_right, sf_right =  subtree(right_edge, prob, model, data, descendants, Ntip, K, elt)
    end

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

function subtree(edge_index, prob, model, data, descendants, Ntip, K, elt)
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
        
        local D_left, D_right
        Threads.@sync begin ## wait for both subtrees to finish before proceeding
            ## spawn a new task for the left subtree, so that it can be computed in parallel
            Threads.@spawn D_left, sf_left = subtree(left_edge, prob, model, data, descendants, Ntip, K, elt)
            D_right, sf_right = subtree(right_edge, prob, model, data, descendants, Ntip, K, elt)
        end

       
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

#=
function postorder_async(model::SSE, data::SSEdata, E)
    ## Pre-compute descendants in hashtable
    descendants = make_descendants(data)

    K = number_of_states(model)

    ## Postorder traversal: computing the branch probabilities through time
    nrows = size(data.edges, 1)
    Ntip = length(data.tiplab)
    
    ## Storing the solution at the end of the branch
    elt = eltype(model)
    D_ends = zeros(elt, nrows, K)

    ## Storing the scaling factors
    sf = zeros(elt, nrows)

    pD = (model.λ, model.μ, model.η, K, E)
    u0 = ones(elt, K)
    tspan = (0.0, 1.0)
    ode = backward_prob(model)
    prob = OrdinaryDiffEq.ODEProblem{true}(ode, u0, tspan, pD)

    root_index = Ntip+1
    left_edge, right_edge = descendants[root_index]

    Threads.@sync begin
        Threads.@spawn subtree!(left_edge, D_ends, sf, prob, model, data, descendants, Ntip, K, elt)
        Threads.@spawn subtree!(right_edge, D_ends, sf, prob, model, data, descendants, Ntip, K, elt)
    end

    return(D_ends, sf)
end

function subtree!(edge_index, D_ends, sf, prob, model, data, descendants, Ntip, K, elt)
    anc, dec = data.edges[edge_index,:]
    alg = OrdinaryDiffEq.Tsit5()  
    node_age = data.node_depth[dec]
    parent_node_age = data.node_depth[anc]
    tspan = (node_age, parent_node_age)  

    if dec <= Ntip
        u0 = ones(elt, K) .* data.sampling_probability
    else
        left_edge, right_edge = descendants[dec]
        
        Threads.@sync begin ## wait for both subtrees to finish before proceeding
            ## spawn a new task for the left subtree, so that it can be computed in parallel
            Threads.@spawn subtree!(left_edge, D_ends, sf, prob, model, data, descendants, Ntip, K, elt)
            subtree!(right_edge, D_ends, sf, prob, model, data, descendants, Ntip, K, elt)
        end

        D_left = D_ends[left_edge,:]
        D_right = D_ends[right_edge,:]
       
        λt = get_speciation_rates(model, node_age) 
        u0 = D_left .* D_right .* λt
    end

    prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
    sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg, save_everystep = false, reltol = 1e-3)
    
    sol = sol.u[end]
    c = sum(sol)
    sol = sol ./ c
    D_ends[edge_index,:] = sol
    if c > 0.0
        sf[edge_index] += log(c)
    end
end

=#
