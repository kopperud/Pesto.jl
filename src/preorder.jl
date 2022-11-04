export preorder

function preorder(model, data, E, D_ends; verbose = false, alg = DifferentialEquations.Tsit5())
    ## Preorder pass, compute `F(t)`
    k = length(model.λ)
    ntips = length(data.tiplab)
    i_not_js = [setdiff(1:k, i) for i in 1:k]
    
    root_node = length(data.tiplab)+1
    left_root_edge = findall(data.edges[:,1] .== root_node)[1]

    nrows = size(data.edges, 1)
    ## Store the numerical solution of F at the end of the branch
    F_ends = zeros(typeof(model.λ[1]), nrows, k)
    ## Store the whole `F(t)` per branch
    Fs = Dict()

    pF = [model.λ, model.μ, model.η, i_not_js, k, E]

    if verbose
        prog = ProgressMeter.Progress(length(data.po), "Preorder pass")
    end

    for i in reverse(data.po)
        parent = data.edges[i,1]
        child = data.edges[i,2]
        
        ## if root
        if parent == root_node
            root_children = findall(data.edges[:,1] .== root_node)
            other_child = setdiff(root_children, i)[1]

            F_start = D_ends[other_child,:] .* model.λ
        else
            parent_edge = findall(data.edges[:,2] .== parent)[1]
            children = findall(data.edges[:,1] .== parent)
            other_child = setdiff(children, i)[1]

            F_start = F_ends[parent_edge,:] .* model.λ .* D_ends[other_child,:]
            F_start = F_start ./ sum(F_start) ## Normalize, because these numbers can get very tiny (1E-10)
        end

        node_age = data.node_depth[child]
        parent_node = parental_node(child, data)
        parent_node_age = data.node_depth[parent_node]
        tspan = (parent_node_age, node_age)

        u0 = F_start
        prob = DifferentialEquations.ODEProblem(forward_prob, u0, tspan, pF)
        #sol = DifferentialEquations.solve(prob, alg, save_everystep = false)[end]
        sol = DifferentialEquations.solve(prob, alg, isoutofdomain = (u,p,t)->any(x->x<0,u))
        Fs[i] = sol
        sol = sol[end]

        F_ends[i,:] = sol

        if verbose
            ProgressMeter.next!(prog)
        end
    end

    return(Fs, F_ends)
end

export ancestral_state_probs
function ancestral_state_probs(data, model, D_ends, F_ends)
    k = length(model.λ)
    root_node = length(data.tiplab)+1
    ntips = length(data.tiplab)

    ASP = zeros(ntips-1, k)

    for node in root_node : maximum(data.edges)
        children = findall(data.edges[:,1] .== node)

        if node == root_node
            asp = D_ends[children[1],:] .* D_ends[children[2],:] .* model.λ
        else
            edge_idx = findall(data.edges[:,2] .== node)[1]
            asp =  D_ends[children[1],:] .* D_ends[children[2],:] .* model.λ .* F_ends[edge_idx,:]
        end
        ## Assign and normalize
        ASP[node - ntips,:] = asp ./ sum(asp)
    end
    return(ASP)
end


