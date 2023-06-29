export postorder_nosave

notneg(u,p,t) = any(x->x<0,u)

"""
    postorder_nosave(model::SSEconstant, data::SSEdata, E, alg = OrdinaryDiffEq.Tsit5())

TBW
"""
function postorder_nosave(model::SSEconstant, data::SSEdata, E, alg = OrdinaryDiffEq.Tsit5())
    ## Pre-compute descendants in hashtable
    descendants = make_descendants(data)
    ancestors = make_ancestors(data)

    n = length(model.λ)

    ## Postorder traversal: computing the branch probabilities through time
    nrows = size(data.edges, 1)
    Ntip = length(data.tiplab)
    
    ## Storing the solution at the end of the branch
    D_ends = zeros(typeof(model.η), nrows, n)
    ## Storing the scaling factors
    sf = zeros(typeof(model.η), nrows)

    pD = (model.λ, model.μ, model.η, n, E)
    D0 = repeat([1.0], n)
    u0 = typeof(model.η).(D0)
    tspan = (0.0, 1.0)
    prob = OrdinaryDiffEq.ODEProblem(backward_prob, u0, tspan, pD)
 
    #Threads.@threads for i in data.po
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

            u0 = typeof(model.λ[1]).(D)

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

    index_partitions = partition_postorder_indices(data)

    for indices in index_partitions[2:end]
        #Threads.@threads for node_idx in indices
        for node_idx in indices
            i = ancestors[node_idx]
            anc, dec = data.edges[i,:]
            if dec > Ntip
                
                left_edge, right_edge = descendants[dec]            
                node_age = data.node_depth[dec]

                D_left = D_ends[left_edge,:]
                D_right = D_ends[right_edge,:]

                D = D_left .* D_right .* model.λ
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
    end
    return(D_ends, sf)
end

