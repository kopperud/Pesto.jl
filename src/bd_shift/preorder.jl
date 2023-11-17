export preorder

function preorder(model::SSE, data::SSEdata, E, Ds; alg = OrdinaryDiffEq.Tsit5())
    ## Precompute ancestor edges
    ancestors = make_ancestors(data)
    descendants = make_descendants(data)

    ## Preorder pass, compute `F(t)`
    k = number_of_states(model)
    elt = eltype(model)
    
    root_node = length(data.tiplab)+1  

    nrows = size(data.edges, 1)
    ## Store the whole `F(t)` per branch
    Fs = Dict()

    pF = (model.λ, model.μ, model.η, k, E)
    ode = forward_prob(model)

    for m in reverse(data.po)
        anc = data.edges[m,1]
        dec = data.edges[m,2]
        
        ## if root
        if anc == root_node
            F_parent = ones(elt, k)
            left_edge, right_edge = descendants[root_node]
            root_age = maximum(data.node_depth)
            λroot = get_speciation_rates(model, root_age)
            D_parent = Ds[left_edge][end] .* Ds[right_edge][end] .* λroot
        else
            parent_edge = ancestors[anc]
            F_parent = Fs[parent_edge][end]
            D_parent = Ds[parent_edge][1]
        end
        Dm = Ds[m][end]

        F_start = D_parent .* F_parent ./ Dm
        F_start = F_start ./ sum(F_start) ## Normalize, because these numbers can get very tiny (1E-10)

        node_age = data.node_depth[dec]
        parent_node = parental_node(dec, data)
        parent_node_age = data.node_depth[parent_node]
        tspan = (parent_node_age, node_age)

        u0 = F_start
        prob = OrdinaryDiffEq.ODEProblem(ode, u0, tspan, pF)
        sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = (u,p,t)->any(x->x<0,u))
        Fs[m] = sol
    end

    return(Fs)
end
