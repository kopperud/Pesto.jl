export postorder

function postorder(model::Model, data::SSEdata)
    ## Pre-compute descendants in hashtable
    descendants = make_descendants(data)
    K = number_of_states(model)

    ## Postorder traversal: computing the branch probabilities through time
    nrows = size(data.edges, 1)
    Ntip = length(data.tiplab)
    elt = eltype(model)
    alg = OrdinaryDiffEq.Tsit5()

    ## Storing the Ds
    us = Dict{Int64, OrdinaryDiffEq.ODESolution}()
    ## Storing the solution at the end of the branch
    u_ends = zeros(elt, nrows, K, 2)

    p = (model, K)
    u0 = ones(elt, K, 2)
    tspan = (0.0, 1.0)
    ode = backward_prob(model)
    prob = OrdinaryDiffEq.ODEProblem(ode, u0, tspan, p)
 
    for m in data.po
        anc, dec = data.edges[m,:]

        if dec < Ntip+1

            species = data.tiplab[dec]
            trait_value = data.trait_data[species]

            if trait_value == "?" ## If we don't know or or didn't observe the trait
                D = repeat([1.0], K) .* data.sampling_probability                
            else ## If we observed the trait and measured it 
                trait_idx = convert.(Float64, trait_value .== data.state_space)
                D = trait_idx .* data.sampling_probability
            end
            E = repeat([1.0 - data.sampling_probability], K)

            u = hcat(E, D)
            u0 = elt.(u)

            node_age = data.node_depth[dec]
            parent_node_age = data.node_depth[anc]
            tspan = (node_age, parent_node_age)

            prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
            sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg)
            us[m] = sol
            E = sol.u[end][:,1]
            D = sol.u[end][:,2]

            c = sum(D)
            D = D ./ c
            u = hcat(E, D)

            u_ends[m,:,:] = u

            logc = log(c)
        end
    end

    for m in data.po

        anc, dec = data.edges[m,:]
        if dec > Ntip

            left_edge, right_edge = descendants[dec]            
            node_age = data.node_depth[dec]

            u_left = u_ends[left_edge,:,:]
            u_right = u_ends[right_edge,:,:]

            D_left = u_left[:,2]
            D_right = u_right[:,2]
            E_left = u_left[:,1]
                       
            λt = get_speciation_rates(model, node_age)
            D = D_left .* D_right .* λt 
            u0 = hcat(E_left, D)

            parent_node_age = data.node_depth[anc]
            tspan = (node_age, parent_node_age)

            prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
            sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg)

            us[m] = sol

            E = sol.u[end][:,1]
            D = sol.u[end][:,2]

            c = sum(D)
            D = D ./ c
            u = hcat(E, D)

            u_ends[m,:,:] = u
        end
    end
    return(us)
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
        #E::OrdinaryDiffEq.ODESolution,
       )
   
    #E = extinction_probability(model, root);

    K = number_of_states(model)
    ode = backward_prob(model)
    p = (model, K)
    tspan = (0.0, 1.0)
    u0 = ones(K,2)

    prob = OrdinaryDiffEq.ODEProblem{true}(ode, u0, tspan, p)

    height = treeheight(root);

    ## find sampling probability
    ## assume it is equal in all the species
    leftmost_tip = find_one_extant_tip(root)
    sampling_probability = leftmost_tip.sampling_probability

    Ds = Dict{Int64, OrdinaryDiffEq.ODESolution}()

    D_root = postorder!(model, root, prob, sampling_probability, height, Ds)

    return(Ds)
end

## internal node
function postorder!(
        model::Model, 
        node::T, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64, 
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
        )  where {T <: BranchingEvent}

    branch_left, branch_right = node.children

    local D_left, D_right

    ## note this is not thread safe
    u_left = postorder!(model, branch_left, prob, sampling_probability, time, Ds)
    u_right = postorder!(model, branch_right, prob, sampling_probability, time, Ds)

    D_left = u_left[:,2]
    D_right = u_right[:,2]
    E_left = u_left[:,1]

    D = D_left .* D_right .* model.λ
    c = sum(D)
    D = D ./ c

    u = hcat(E_left, D)

    return(u)
end

## along a branch
function postorder!(
        model::Model, 
        branch::Branch, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64,
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
    )
    child_node = branch.outbounds
    t_old = time 
    t_young = time - branch.time

    u0 = postorder!(model, child_node, prob, sampling_probability, t_young, Ds)

    tspan = (t_young, t_old)
    prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5(), isoutofdomain = notneg, save_everystep = true, reltol = 1e-3)

    u = sol.u[end]

    c = sum(u[:,2])
    u[:,2] .= u[:,2] ./ c

    Ds[branch.index] = sol

    return(u)
end


## for a tip
function postorder!(
        model::Model, 
        tip::ExtantTip, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64,
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
    )
    elt = eltype(model)
    K = number_of_states(model)

    E = ones(elt, K) .- tip.sampling_probability
    D = zeros(elt, K) .+ tip.sampling_probability
    #sf = 0.0

    u = hcat(E, D)

    return(u)
end


## for a tip
function postorder!(
        model::Model, 
        tip::FossilTip, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64,
        Ds::Dict{Int64, OrdinaryDiffEq.ODESolution},
    )
    elt = eltype(model)
    K = number_of_states(model)

    r = 0.0 ## fixed

    ψ = get_fossilization_rate(model, time)

    Et = extinction_probability(model, sampling_probability, time)

    D = ψ .* Et .* (1 - r)
    D += (r .* ψ)

    u = hcat(Et, D)

    #sf = 0.0

    return(u)
end




