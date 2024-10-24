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

model = BDSconstant(λ, μ, η)

E = extinction_probability(model, data);
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
function postorder_async(model::Model, data::SSEdata, E)
    ## Pre-compute descendants in hashtable
    descendants = make_descendants(data)

    K = number_of_states(model)

    ## Postorder traversal: computing the branch probabilities through time
    nrows = size(data.edges, 1)
    Ntip = length(data.tiplab)
    
    ## Storing the solution at the end of the branch
    elt = eltype(model)

    ## Storing the scaling factors
    pD = (model, K, E)

    u0 = ones(elt, K)
    tspan = (0.0, 1.0)
    ode = backward_prob(model)
    prob = OrdinaryDiffEq.ODEProblem{true}(ode, u0, tspan, pD)

    root_index = Ntip+1
    root_age = data.node_depth[root_index]

    left_edge, right_edge = descendants[root_index]

    local sf_left, sf_right
    local D_left, D_right
    Threads.@sync begin
        Threads.@spawn D_left, sf_left = subtree(left_edge, prob, model, data, descendants, Ntip, K, elt)
        D_right, sf_right =  subtree(right_edge, prob, model, data, descendants, Ntip, K, elt)
    end

    λroot = get_speciation_rates(model, root_age)
    D = D_left .* D_right .* λroot
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

######################################################
##
##   compute it with the other tree format
##
######################################################

## the root
function postorder_async(
        model::Model,
        root::Root,
        E::OrdinaryDiffEq.ODESolution,
       )
   
    E = extinction_probability(model, root);

    elt = eltype(model)
    K = number_of_states(model)
    ode = backward_prob(model)
    pD = (model, K, E)
    tspan = (0.0, 1.0)
    u0 = ones(elt, K)

    prob = OrdinaryDiffEq.ODEProblem{true}(ode, u0, tspan, pD)

    height = treeheight(root);

    D, sf = postorder_async(model, root, prob, height, E)
    return(D, sf)
end

## internal node
function postorder_async(
        model::Model, 
        node::T, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64, 
        E,
        )  where {T <: BranchingEvent}

    branch_left, branch_right = node.children

    local D_left, D_right
    local sf_left, sf_right

    Threads.@sync begin
        Threads.@spawn D_left, sf_left = postorder_async(model, branch_left, prob, time, E)
        D_right, sf_right = postorder_async(model, branch_right, prob, time, E)
    end

    D = D_left .* D_right .* model.λ
    c = sum(D)
    sf = sf_left + sf_right + log(c)
    D = D ./ c

    return(D, sf)
end

## along a branch
function postorder_async(
        model::Model, 
        branch::Branch, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64,
        E,
    )
    child_node = branch.outbounds
    t_old = time 
    t_young = time - branch.time

    D0, sf = postorder_async(model, child_node, prob, t_young, E)

    tspan = (t_young, t_old)
    prob = OrdinaryDiffEq.remake(prob, u0 = D0, tspan = tspan)
    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5(), isoutofdomain = notneg, 
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


    return(D, sf)
end


## for a tip
function postorder_async(
        model::Model, 
        tip::ExtantTip, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64,
        E,
    )
    elt = eltype(model)
    K = number_of_states(model)

    D = ones(elt, K) .* tip.sampling_probability
    sf = 0.0

    return(D, sf)
end

## for a tip
function postorder_async(
        model::Model, 
        tip::FossilTip, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64,
        E,
    )
    elt = eltype(model)
    K = number_of_states(model)

    ## probability that sampled lineages immediately go extinct
    r = 0.0

    ψ = get_fossilization_rate(model, time)
    Et = E(time)

    ## MacPherson et al. (2022) Sys Bio
    D = r .* ψ .+ (1.0 - r) .* ψ .* Et
    sf = 0.0

    return(D, sf)
end




