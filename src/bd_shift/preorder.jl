export preorder

function preorder(
        model::Model, 
        data::SSEdata, 
        post::Dict{Int64, OrdinaryDiffEq.ODESolution}; 
        alg = OrdinaryDiffEq.Tsit5(),
        condition = [:survival, :mrca],
    )
    ## Precompute ancestor edges
    ancestors = make_ancestors(data)
    descendants = make_descendants(data)

    ## Preorder pass, compute `F(t)`
    K = number_of_states(model)
    elt = eltype(model)
    
    root_node = length(data.tiplab)+1  
    nrows = size(data.edges, 1)

    ## Store the whole `F(t)` per branch
    pre = Dict{Int64, OrdinaryDiffEq.ODESolution}()

    p = (model, K)
    ode = forward_prob(model)
    tspan = [0.0, 1.0]
    u0 = zeros(Float64, K, 2) 
    prob = OrdinaryDiffEq.ODEProblem(ode, u0, tspan, p)

    for m in reverse(data.po)
        anc = data.edges[m,1]
        dec = data.edges[m,2]
       
        ## if root
        if anc == root_node
            ## initialize F as the prior probability
            ## of the base distribution
            F_parent = ones(elt, K) ./ K

            left_edge, right_edge = descendants[root_node]
            root_age = maximum(data.node_depth)
            λroot = get_speciation_rates(model, root_age)

            ## potentially impose conditions
            if :survival in condition
                E = post[left_edge].u[end][:,1]
                nonextinct = (1.0 .- E).^2
                F_parent = F_parent ./ nonextinct
            end

            if :mrca in condition
                λroot = get_speciation_rates(model, root_age)
                F_parent = F_parent ./ λroot
            end

            D_parent = post[left_edge].u[end][:,2] .* post[right_edge].u[end][:,2] .* λroot
            E_parent = post[left_edge].u[end][:,1]
        else
            parent_edge = ancestors[anc]
            F_parent = pre[parent_edge].u[end][:,2]
            D_parent = post[parent_edge].u[1][:,2]
            E_parent = post[parent_edge].u[1][:,1]
        end
        ## D(t) on this branch
        Dm = post[m].u[end][:,2]

        F_start = D_parent .* F_parent ./ Dm
        ## Normalize for numerical stability 
        F_start = F_start ./ sum(F_start) 


        parent_node = parental_node(dec, data) 

        ## node age (youngest) 
        node_age = data.node_depth[dec] 
        ## parent node age (oldest)
        parent_node_age = data.node_depth[parent_node] 
        tspan = (parent_node_age, node_age)

        u0 = hcat(E_parent, F_start)
        prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
        sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg)
        pre[m] = sol
    end

    return(pre)
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
        post::Dict{Int64, OrdinaryDiffEq.ODESolution};
        condition = [:mrca, :survival],
       )
   
    #E = extinction_probability(model, root);

    elt = eltype(model)
    K = number_of_states(model)
    ode = forward_prob(model)
    p = (model, K)
    tspan = (0.0, 1.0)
    u0 = ones(elt, K, 2)

    prob = OrdinaryDiffEq.ODEProblem{true}(ode, u0, tspan, p)

    height = treeheight(root);

    pre = Dict{Int64, OrdinaryDiffEq.ODESolution}()

    left_index = root.children[1].index
    right_index = root.children[2].index

    F_root = ones(K) ./ K
    D_left = post[left_index].u[end][:,2] 
    D_right = post[right_index].u[end][:,2] 
    λroot = get_speciation_rates(model, height)
    D_root = D_left .* D_right .* λroot
    S_root = F_root .* D_root

    ## potentially impose conditions
    if :survival in condition
        E = post[left_index].u[end][:,1]
        nonextinct = (1.0 .- E).^2
        S_root = S_root ./ nonextinct
    end

    if :mrca in condition
        root_age = treeheight(root)
        λroot = get_speciation_rates(model, root_age)
        S_root = S_root ./ λroot
    end

    #normalize
    S_root = S_root ./ sum(S_root)

    preorder!(model, root, prob, height, post, pre, S_root)

    return(pre)
end

## internal node
function preorder!(
        model::Model, 
        node::T, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64, 
        post::Dict{Int64, OrdinaryDiffEq.ODESolution},
        pre::Dict{Int64, OrdinaryDiffEq.ODESolution},
        S_node::Vector{Float64}
        )  where {T <: BranchingEvent}

    left_branch_index = node.children[1].index
    right_branch_index = node.children[2].index

    branch_left, branch_right = node.children

    # not thread safe because Dict is not thread safe
    preorder!(model, branch_left, prob, time, post, pre, S_node)
    preorder!(model, branch_right, prob, time, post, pre, S_node)
end

## sampled ancestor
function preorder!(
        model::Model, 
        node::SampledAncestor, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64, 
        post::Dict{Int64, OrdinaryDiffEq.ODESolution},
        pre::Dict{Int64, OrdinaryDiffEq.ODESolution},
        S_node::Vector{Float64}
        )

    child = node.child
    S_node = S_node ./ model.ψ
    S_node = S_node ./ sum(S_node)

    preorder!(model, child, prob, time, post, pre, S_node)
end


## along a branch
function preorder!(
        model::Model, 
        branch::Branch, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64,
        post::Dict{Int64, OrdinaryDiffEq.ODESolution},
        pre::Dict{Int64, OrdinaryDiffEq.ODESolution},
        S_parent::Vector{Float64},
    )
    child_node = branch.outbounds
    t_old = time 
    t_young = time - branch.time

   
    D0 = post[branch.index].u[end][:,2]
    F0 = S_parent ./ D0
    F0 = F0 ./ sum(F0)
    E0 = post[branch.index].u[end][:,1]

    u0 = hcat(E0, F0)

    tspan = (t_old, t_young)
    prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)

    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5(), isoutofdomain = notneg, save_everystep = true, reltol = 1e-3)

    push!(pre, branch.index => sol)

    u = sol.u[end]
    F_young = u[:,2]
    D_young = post[branch.index].u[1][:,2]

    S_young = F_young .* D_young
    c = sum(S_young) 
    S_young = S_young ./ c

    preorder!(model, child_node, prob, t_young, post, pre, S_young)
end

## for a tip
function preorder!(
        model::Model, 
        tip::T, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64,
        post::Dict{Int64, OrdinaryDiffEq.ODESolution},
        pre::Dict{Int64, OrdinaryDiffEq.ODESolution},
        S_parent::Vector{Float64},
    ) where {T <: AbstractTip}
    nothing
end

