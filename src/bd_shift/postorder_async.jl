#asd
export postorder_async
## async log likelihood calculation

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
    prob = OrdinaryDiffEq.ODEProblem(ode, u0, tspan, pD)

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
        u0 = ones(elt, K) .* data.ρ
    else
        left_edge, right_edge = descendants[dec]
        
        Threads.@sync begin ## wait for both subtrees to finish before proceeding
            ## spawn a new task for both subtrees, so that they can be computed in parallel
            Threads.@spawn subtree!(left_edge, D_ends, sf, prob, model, data, descendants, Ntip, K, elt)
            Threads.@spawn subtree!(right_edge, D_ends, sf, prob, model, data, descendants, Ntip, K, elt)
        end

        D_left = D_ends[left_edge,:]
        D_right = D_ends[right_edge,:]
       
        λt = get_speciation_rates(model, node_age) 
        u0 = D_left .* D_right .* λt
    end

    prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
    sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg, save_everystep = false)
    
    sol = sol[end]
    c = sum(sol)
    sol = sol ./ c
    D_ends[edge_index,:] = sol
    if c > 0.0
        sf[edge_index] += log(c)
    end
end
