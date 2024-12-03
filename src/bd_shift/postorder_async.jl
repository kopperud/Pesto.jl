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
function postorder_async(model::Model, data::SSEdata)
    ## Pre-compute descendants in hashtable
    descendants = make_descendants(data)

    K = number_of_states(model)

    ## Postorder traversal: computing the branch probabilities through time
    nrows = size(data.edges, 1)
    Ntip = length(data.tiplab)
    
    ## Storing the solution at the end of the branch
    elt = eltype(model)

    ## Storing the scaling factors
    p = (model, K)

    u0 = ones(elt, K, 2)
    tspan = (0.0, 1.0)
    ode = backward_prob(model)
    prob = OrdinaryDiffEq.ODEProblem{true}(ode, u0, tspan, p)

    root_index = Ntip+1
    root_age = data.node_depth[root_index]

    left_edge, right_edge = descendants[root_index]

    local sf_left, sf_right
    local u_left, u_right
    Threads.@sync begin
        Threads.@spawn u_left, sf_left = subtree(left_edge, prob, model, data, descendants, Ntip, K, elt)
        u_right, sf_right =  subtree(right_edge, prob, model, data, descendants, Ntip, K, elt)
    end

    D_left = u_left[:,2]
    D_right = u_right[:,2]
    E_left = u_left[:,1]

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

    u = hcat(E_left, D)

    return(u, sf)
end

function subtree(edge_index, prob, model, data, descendants, Ntip, K, elt)
    anc, dec = data.edges[edge_index,:]
    alg = OrdinaryDiffEq.Tsit5()  
    node_age = data.node_depth[dec]
    parent_node_age = data.node_depth[anc]
    tspan = (node_age, parent_node_age)  

    local sf_left, sf_right
    if dec <= Ntip
        u0 = zeros(elt, K, 2) 
        u0[:,1] .= 1 - data.sampling_probability
        u0[:,2] .= data.sampling_probability

        sf_left = 0.0
        sf_right = 0.0
    else
        left_edge, right_edge = descendants[dec]
        
        local u_left, u_right
        Threads.@sync begin ## wait for both subtrees to finish before proceeding
            ## spawn a new task for the left subtree, so that it can be computed in parallel
            Threads.@spawn u_left, sf_left = subtree(left_edge, prob, model, data, descendants, Ntip, K, elt)
            u_right, sf_right = subtree(right_edge, prob, model, data, descendants, Ntip, K, elt)
        end

        D_left = u_left[:,2]
        D_right = u_right[:,2]
        E_left = u_left[:,1]
         
        λt = get_speciation_rates(model, node_age) 
        D0 = D_left .* D_right .* λt
        E0 = E_left
        u0 = hcat(E0, D0)
    end

    sf = sf_left + sf_right
    if isinf(sf) | isnan(sf) | in(u0, NaN)
        u = ones(elt, K, 2)
        sf -= Inf
    else
        prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
        sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg, 
                                   save_everystep = false, reltol = 1e-3)
        u = sol.u[end]
        c = sum(u[:,2])
        u[:,2] = u[:,2] ./ c

        if c > 0.0
            sf += log(c)
        else
            sf -= Inf
        end

        if !(sol.retcode == OrdinaryDiffEq.ReturnCode.Success)
            sf -= Inf
        end
    end

    return(u, sf)
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
    )
   
    #E = extinction_probability(model, root);

    elt = eltype(model)
    K = number_of_states(model)
    ode = backward_prob(model)
    p = (model, K)
    tspan = (0.0, 1.0)
    u0 = ones(elt, K, 2)

    prob = OrdinaryDiffEq.ODEProblem{true,SciMLBase.FullSpecialize}(ode, u0, tspan, p)

    height = treeheight(root);

    ## find sampling probability
    ## assume it is equal in all the species
    leftmost_tip = find_one_extant_tip(root)
    sampling_probability = leftmost_tip.sampling_probability

    u, sf = postorder_async(model, root, prob, sampling_probability, height)
    return(u, sf)
end

## internal node
function postorder_async(
        model::Model, 
        node::N, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64, 
        ) where {N <: BranchingEvent}

    branch_left, branch_right = node.children

    local u_left, u_right
    local sf_left, sf_right

    Threads.@sync begin
        Threads.@spawn u_left, sf_left = postorder_async(model, branch_left, prob,sampling_probability,  time)
        u_right, sf_right = postorder_async(model, branch_right, prob,sampling_probability,  time)
    end

    D_left = u_left[:,2]
    D_right = u_right[:,2]
    E_left = u_left[:,1]

    D = D_left .* D_right .* model.λ
    c = sum(D)
    sf = sf_left + sf_right + log(c)
    D = D ./ c

    u = hcat(E_left, D)

    return(u, sf)
end


## along a branch
function postorder_async(
        model::Model, 
        branch::Branch, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64,
    )
    child_node = branch.outbounds

    t_old = time 
    t_young = time - branch.time

    u0, sf = postorder_async(model, child_node, prob, sampling_probability, t_young)

    tspan = (t_young, t_old)
    prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5(), isoutofdomain = notneg, save_everystep = false, reltol = 1e-3)

    u = sol.u[end]
    c = sum(u, dims = 1)[2]
    u[:,2] = u[:,2] ./ c

    if c > 0.0
        sf += log(c)
    else
        sf -= Inf
    end

    if !(sol.retcode == OrdinaryDiffEq.ReturnCode.Success)
        sf -= Inf
    end


    return(u, sf)
end


## for a tip
function postorder_async(
        model::Model, 
        tip::ExtantTip, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64,
    )

    @assert abs(time .- 0) < 0.001
    elt = eltype(model)
    K = number_of_states(model)

    E = ones(elt, K) .- tip.sampling_probability
    D = zeros(elt, K) .+ tip.sampling_probability
    sf = 0.0

    u = hcat(E, D)

    return(u, sf)
end

## for a tip
function postorder_async(
        model::Model, 
        tip::FossilTip, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64,
    )

    elt = eltype(model)
    K = number_of_states(model)

    ## probability that sampled lineages immediately go extinct
    r = 0.0

    ψ = get_fossilization_rate(model, time)
    ## Problem:
    ## How do we know the extinction probability E(t) for this 
    ## fossil tip, if we don't know the extant sampling fraction, 
    ## or if we assume that the extant sampling probability 
    ## is unequal across species?
    Et = extinction_probability(model, sampling_probability, time)

    ## MacPherson et al. (2022) Sys Bio
    D = r .* ψ .+ (1.0 - r) .* ψ .* Et
    sf = 0.0

    u = hcat(Et, D)

    return(u, sf)
end




