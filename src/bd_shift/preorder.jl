export preorder

function preorder(model::Model, data::SSEdata, E, Ds; alg = OrdinaryDiffEq.Tsit5())
    ## Precompute ancestor edges
    ancestors = make_ancestors(data)
    descendants = make_descendants(data)

    ## Preorder pass, compute `F(t)`
    K = number_of_states(model)
    elt = eltype(model)
    
    root_node = length(data.tiplab)+1  

    nrows = size(data.edges, 1)
    ## Store the whole `F(t)` per branch
    Fs = Dict()

    #pF = (model.λ, model.μ, model.η, K, E)
    pF = (model, K, E)
    ode = forward_prob(model)
    tspan = [0.0, 1.0]
    u0 = [1.0, 1.0]
    prob = OrdinaryDiffEq.ODEProblem(ode, u0, tspan, pF)

    for m in reverse(data.po)
        anc = data.edges[m,1]
        dec = data.edges[m,2]
       
        ## if root
        if anc == root_node
            F_parent = ones(elt, K) ./ K

            left_edge, right_edge = descendants[root_node]
            root_age = maximum(data.node_depth)
            λroot = get_speciation_rates(model, root_age)
            D_parent = Ds[left_edge].u[end] .* Ds[right_edge].u[end] .* λroot
        else
            parent_edge = ancestors[anc]
            F_parent = Fs[parent_edge].u[end]
            D_parent = Ds[parent_edge].u[1]
        end
        Dm = Ds[m].u[end]

        F_start = D_parent .* F_parent ./ Dm
        F_start = F_start ./ sum(F_start) ## Normalize, because these numbers can get very tiny (1E-10)

        parent_node = parental_node(dec, data) 

        node_age = data.node_depth[dec] ## node age (youngest) 
        parent_node_age = data.node_depth[parent_node] ## parent node age (oldest)
        tspan = (parent_node_age, node_age)

        u0 = F_start
        prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
        sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg)
        Fs[m] = sol
    end

    return(Fs)
end


######################################################
##
##   compute it with the other tree format
##
######################################################

## the root
function preorder(
        model::Model,
        root::Root,
        E::OrdinaryDiffEq.ODESolution,
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
       )
   
    E = extinction_probability(model, root);

    elt = eltype(model)
    K = number_of_states(model)
    ode = forward_prob(model)
    pF = (model, K, E)
    tspan = (0.0, 1.0)
    u0 = ones(elt, K)

    prob = OrdinaryDiffEq.ODEProblem{true}(ode, u0, tspan, pF)

    height = treeheight(root);

    ## 

    Fs = Dict{Int64, OrdinaryDiffEq.ODESolution}()

    F_root = ones(K) 
    preorder!(model, root, prob, height, E, Ds, Fs, F_root)

    return(Fs)
end

## internal node
function preorder!(
        model::Model, 
        node::T, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64, 
        E::OrdinaryDiffEq.ODESolution,
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
        Fs::Dict{Int64, OrdinaryDiffEq.ODESolution},
        F_node::Vector{Float64}
        )  where {T <: BranchingEvent}

    # calculate D(t) at this node (or the parent edge of this node)
    left_branch_index = node.children[1].index
    right_branch_index = node.children[2].index

    D_left = Ds[left_branch_index](time)
    D_right = Ds[right_branch_index](time)
    λt = get_speciation_rates(model, time)
    D_node = D_left .* D_right .* λt

    # calculate S(t)
    S_node = F_node .* D_node
    S_node = S_node ./ sum(S_node) ## normalize

    branch_left, branch_right = node.children

    # not thread safe because Dict is not thread safe
    preorder!(model, branch_left, prob, time, E, Ds, Fs, S_node)
    preorder!(model, branch_right, prob, time, E, Ds, Fs, S_node)
end

## along a branch
function preorder!(
        model::Model, 
        branch::Branch, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64,
        E::OrdinaryDiffEq.ODESolution,
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
        Fs::Dict{Int64, OrdinaryDiffEq.ODESolution},
        S_parent::Vector{Float64},
    )
    child_node = branch.outbounds
    t_old = time 
    t_young = time - branch.time

    
    F0 = S_parent ./ Ds[branch.index](t_old)
    F0 = F0 ./ sum(F0)

    tspan = (t_old, t_young)
    prob = OrdinaryDiffEq.remake(prob, u0 = F0, tspan = tspan)

    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5(), isoutofdomain = notneg, 
                               save_everystep = true, reltol = 1e-3)
    #sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5())

    #Fs[branch.index] = sol
    push!(Fs, branch.index => sol)

    F_young = sol.u[end]
    c = sum(F_young) 
    F_young = F_young ./ c

    preorder!(model, child_node, prob, t_young, E, Ds, Fs, F_young)
end

## for a tip
function preorder!(
        model::Model, 
        tip::T, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64,
        E::OrdinaryDiffEq.ODESolution,
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
        Fs::Dict{Int64, OrdinaryDiffEq.ODESolution},
        S_parent::Vector{Float64},
    ) where {T <: AbstractTip}
    nothing
end



