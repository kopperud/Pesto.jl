export postorder

function postorder(model::Model, data::SSEdata, E; alg = OrdinaryDiffEq.Tsit5())
    ## Pre-compute descendants in hashtable
    descendants = make_descendants(data)

    n = number_of_states(model)

    ## Postorder traversal: computing the branch probabilities through time
    nrows = size(data.edges, 1)
    Ntip = length(data.tiplab)
    elt = eltype(model)

    ## Storing the Ds
    Ds = Dict()
    ## Storing the solution at the end of the branch
    D_ends = zeros(elt, nrows, n)
#    D_begins = zeros(elt, nrows, n)
    ## Storing the scaling factors
    sf = zeros(elt, nrows)

    #pD = (model.λ, model.μ, model.η, n, E)
    pD = (model, n, E)
    u0 = ones(elt, n)
    tspan = (0.0, 1.0)
    ode = backward_prob(model)
    prob = OrdinaryDiffEq.ODEProblem(ode, u0, tspan, pD)
 
    for m in data.po
        anc, dec = data.edges[m,:]

        if dec < Ntip+1

            species = data.tiplab[dec]
            trait_value = data.trait_data[species]

            if trait_value == "?" ## If we don't know or or didn't observe the trait
                D = repeat([1.0], n) .* data.sampling_probability                
            else ## If we observed the trait and measured it 
                trait_idx = convert.(Float64, trait_value .== data.state_space)
                D = trait_idx .* data.sampling_probability
            end

            u0 = elt.(D)
            #D_begins[i,:] = u0

            node_age = data.node_depth[dec]
            parent_node_age = data.node_depth[anc]
            tspan = (node_age, parent_node_age)

            prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
            sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg)
            Ds[m] = sol
            sol = sol.u[end]

            k = sum(sol)
            sol = sol ./ k
            D_ends[m,:] = sol
            logk = log(k)
            sf[m] = logk
        end
    end

    for m in data.po

        anc, dec = data.edges[m,:]
        if dec > Ntip

            left_edge, right_edge = descendants[dec]            
            node_age = data.node_depth[dec]

            D_left = D_ends[left_edge,:]
            D_right = D_ends[right_edge,:]
            
            λt = get_speciation_rates(model, node_age)
            D = D_left .* D_right .* λt 
            u0 = D

            parent_node_age = data.node_depth[anc]
            tspan = (node_age, parent_node_age)

            prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
            sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg)
            Ds[m] = sol
            sol = sol.u[end]
            k = sum(sol)
            sol = sol ./ k
            D_ends[m,:] = sol
            if k > 0.0
                sf[m] += log(k)
            end
        end
    end
    return(Ds, sf)
end


######################################################
##
##   compute it with the other tree format
##
######################################################

## the root
function postorder(
        model::Model,
        root::Root,
        E::OrdinaryDiffEq.ODESolution,
       )
   
    E = extinction_probability(model, root);

    K = number_of_states(model)
    ode = backward_prob(model)
    pD = (model, K, E)
    tspan = (0.0, 1.0)
    u0 = ones(K)

    prob = OrdinaryDiffEq.ODEProblem{true}(ode, u0, tspan, pD)

    height = treeheight(root);

    Ds = Dict{Int64, OrdinaryDiffEq.ODESolution}()

    D_root = postorder!(model, root, prob, height, E, Ds)

    return(Ds)
end

## internal node
function postorder!(
        model::Model, 
        node::T, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64, 
        E::OrdinaryDiffEq.ODESolution,
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
        )  where {T <: BranchingEvent}

    branch_left, branch_right = node.children

    local D_left, D_right

    ## note this is not thread safe
    D_left = postorder!(model, branch_left, prob, time, E, Ds)
    D_right = postorder!(model, branch_right, prob, time, E, Ds)

    D = D_left .* D_right .* model.λ
    c = sum(D)
    D = D ./ c

    return(D)
end

## along a branch
function postorder!(
        model::Model, 
        branch::Branch, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64,
        E::OrdinaryDiffEq.ODESolution,
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
    )
    child_node = branch.outbounds
    t_old = time 
    t_young = time - branch.time

    D0 = postorder!(model, child_node, prob, t_young, E, Ds)

    tspan = (t_young, t_old)
    prob = OrdinaryDiffEq.remake(prob, u0 = D0, tspan = tspan)
    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5(), isoutofdomain = notneg, 
                                   save_everystep = true, reltol = 1e-3)

    D = sol.u[end]
    c = sum(D) 
    D = D ./ c

    Ds[branch.index] = sol

    return(D)
end


## for a tip
function postorder!(
        model::Model, 
        tip::ExtantTip, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64,
        E::OrdinaryDiffEq.ODESolution,
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
    )
    elt = eltype(model)
    K = number_of_states(model)

    D = ones(elt, K) .* tip.sampling_probability
    sf = 0.0

    return(D)
end


## for a tip
function postorder!(
        model::Model, 
        tip::FossilTip, 
        prob::OrdinaryDiffEq.ODEProblem,
        time::Float64,
        E::OrdinaryDiffEq.ODESolution,
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
    )
    elt = eltype(model)
    K = number_of_states(model)

    r = 0.0 ## fixed


    ψ = get_fossilization_rate(model, time)
    Et = E(time)

    D = ψ .* Et .* (1 - r)
    D += (r .* ψ)

    sf = 0.0

    return(D)
end




