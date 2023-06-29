export preorder

function preorder(model, data, E, D_ends; verbose = false, alg = OrdinaryDiffEq.Tsit5())
    ## Precompute ancestor edges
    ancestors = make_ancestors(data)
    descendants = make_descendants(data)

    ## Preorder pass, compute `F(t)`
    k = length(model.λ)
    
    root_node = length(data.tiplab)+1  

    nrows = size(data.edges, 1)
    ## Store the numerical solution of F at the end of the branch
    F_ends = zeros(typeof(model.λ[1]), nrows, k)
    ## Store the whole `F(t)` per branch
    Fs = Dict()

    pF = (model.λ, model.μ, model.η, k, E)

    if verbose
        prog = ProgressMeter.Progress(length(data.po), "Preorder pass ")
    end

    for i in reverse(data.po)
        anc = data.edges[i,1]
        dec = data.edges[i,2]
        
        ## if root
        if anc == root_node
            root_children = descendants[root_node]           
            other_child = setdiff(root_children, i)[1]

            F_start = D_ends[other_child,:] .* model.λ
        else
            parent_edge = ancestors[anc]
            children = descendants[anc]            
            other_child = setdiff(children, i)[1]

            F_start = F_ends[parent_edge,:] .* model.λ .* D_ends[other_child,:]
            F_start = F_start ./ sum(F_start) ## Normalize, because these numbers can get very tiny (1E-10)
        end

        node_age = data.node_depth[dec]
        parent_node = parental_node(dec, data)
        parent_node_age = data.node_depth[parent_node]
        tspan = (parent_node_age, node_age)

        u0 = F_start
        prob = OrdinaryDiffEq.ODEProblem(forward_prob, u0, tspan, pF)
        sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = (u,p,t)->any(x->x<0,u))
        Fs[i] = sol
        sol = sol[end]

        F_ends[i,:] = sol

        if verbose
            ProgressMeter.next!(prog)
        end
    end

    return(Fs, F_ends)
end