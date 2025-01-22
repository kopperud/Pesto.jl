export postorder_async

## async log likelihood calculation
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

@doc raw"""
    postorder_async(model, tree)

This function does the postorder pass, and returns a tuple of (u, sf).
The matrix `u` is of size K x 2, where K are the number of rate 
categories. 

The first column in `u` represents the extinction probability 
    P(Ψ goes extinct|Z == j, θ),
conditional on that the rate category was j at the root, and
conditional on the parameters of the model θ (i.e. all of the diversification
rate categories λj, μj, ψj, as well as the sampling probability ρ).

The second column `Dj` represents 
    Dj = P(Ψ|Z == j,θ),
i.e. the probability density of observing the reconstructed tree Ψ, conditional
on the same as above.

The scaling factor `sf` is there because we normalize the probability 
densities Dj such that they sum to one.

Example:

```julia
λ = [0.1, 0.2, 0.3]
μ = [0.05, 0.05, 0.05]
η = 0.01

model = BhDhModel(λ, μ, η)

phy = readtree(Pesto.path("bears.tre"))
sampling_probability = 1.0
tree = construct_tree(phy, sampling_probability)

Ds, sf = postorder_async(model, tree)
```
"""
function postorder_async(
        model::Model,
        root::Root,
    )
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

    ψ = get_fossilization_rates(model, time)
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




