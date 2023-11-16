export postorder_nosave

notneg(u,p,t) = any(x->x<0,u)

function eltype(model::SSEconstant)
    return(typeof(model.η))
end
function eltype(model::SSEtimevarying)
    return(typeof(model.η(0.0)))
end

"""
    postorder_nosave(model::SSEconstant, data::SSEdata, E, alg = OrdinaryDiffEq.Tsit5())

TBW
"""
function postorder_nosave(model::SSE, data::SSEdata, E, alg = OrdinaryDiffEq.Tsit5())
    ## Pre-compute descendants in hashtable
    descendants = make_descendants(data)
    ancestors = make_ancestors(data)

    #n = length(model.λ)
    n = number_of_states(model)

    ## Postorder traversal: computing the branch probabilities through time
    nrows = size(data.edges, 1)
    Ntip = length(data.tiplab)
    
    ## Storing the solution at the end of the branch
    elt = eltype(model)
    D_ends = zeros(elt, nrows, n)
    ## Storing the scaling factors
    sf = zeros(elt, nrows)

    pD = (model.λ, model.μ, model.η, n, E)
    #D0 = repeat([1.0], n)
    #u0 = typeof(model.η).(D0)
    u0 = ones(elt, n)
    tspan = (0.0, 1.0)
    ode = backward_prob(model)
    prob = OrdinaryDiffEq.ODEProblem(ode, u0, tspan, pD)
 
    for i in data.po
        anc, dec = data.edges[i,:]
        if dec < Ntip+1
            species = data.tiplab[dec]
            trait_value = data.trait_data[species]

            if trait_value == "?" ## If we don't know or or didn't observe the trait
                D = repeat([1.0], n) .* data.ρ                
            else ## If we observed the trait and measured it 
                trait_idx = convert.(Float64, trait_value .== data.state_space)
                D = trait_idx .* data.ρ
            end

#            u0 = typeof(model.λ[1]).(D)
            u0 = elt.(D)             

            node_age = data.node_depth[dec]
            parent_node_age = data.node_depth[anc]
            tspan = (node_age, parent_node_age)

            prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
            sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg, save_everystep = false)
            sol = sol[end]

            k = sum(sol)
            sol = sol ./ k
            D_ends[i,:] = sol
            logk = log(k)
            sf[i] = logk
        end
    end

    for i in data.po
        anc, dec = data.edges[i,:]
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
            sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = (u,p,t)->any(x->x<0,u), save_everystep = false)
            
            sol = sol[end]
            k = sum(sol)
            sol = sol ./ k
            D_ends[i,:] = sol
            if k > 0.0
                sf[i] += log(k)
            end
        end
    end
    return(D_ends, sf)
end

