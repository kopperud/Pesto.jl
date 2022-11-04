export postorder

function postorder(model::SSEconstant, data::SSEdata; verbose = false, alg = DifferentialEquations.Tsit5())
    ## Compute the extinction probability through time
    n = length(model.λ)
    i_not_js = [setdiff(1:n, i) for i in 1:n]
    pE = [model.λ, model.μ, model.η, i_not_js, n]

    tree_height = maximum(data.node_depth)
    tspan = (0.0, tree_height)

    E0 = repeat([1.0 - data.ρ], n)
    pr = DifferentialEquations.ODEProblem(extinction_prob, E0, tspan, pE);
    E = DifferentialEquations.solve(pr, alg);

    ## Postorder traversal: computing the branch probabilities through time
    nrows = size(data.edges, 1)
    Ntip = length(data.tiplab)

    ## Storing the Ds
    Ds = Dict()
    ## Storing the solution at the end of the branch
    D_ends = zeros(typeof(model.λ[1]), nrows, n)
    ## Storing the scaling factors
    sf = zeros(typeof(model.λ[1]), nrows)

    pD = [model.λ, model.μ, model.η, i_not_js, n, E]
    D0 = repeat([1.0], n)
    u0 = typeof(model.λ[1]).(D0)
    tspan = (0.0, 1.0)
    prob = DifferentialEquations.ODEProblem(backward_prob, u0, tspan, pD)
 
    #for i in 1:nrows
    if verbose
        prog = ProgressMeter.Progress(length(data.po), "Postorder pass")
    end

    for i in data.po
        anc, dec = data.edges[i,:]

        #if is tip
        if dec <= Ntip
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
            parent_node = parental_node(dec, data)
            parent_node_age = data.node_depth[parent_node]
            tspan = (node_age, parent_node_age)
#            times = range(node_age, parent_node_age, length = 100)


            #prob = remake(pr, u0 = u0, tspan = tspan)
#            prob = DifferentialEquations.ODEProblem(backward_prob, u0, tspan, pD)
            prob = DifferentialEquations.remake(prob, u0 = u0, tspan = tspan)
            sol = DifferentialEquations.solve(prob, alg, isoutofdomain = (u,p,t)->any(x->x<0,u))
            Ds[i] = sol
            sol = sol[end]

            k = sum(sol)
            sol = sol ./ k
            D_ends[i,:] = sol
            logk = log(k)
            sf[i] = logk
        else
            left_node, right_node = descendant_nodes(dec, data)
            left_edge = findall(data.edges[:,2] .== left_node)[1]
            right_edge = findall(data.edges[:,2] .== right_node)[1]
            node_age = data.node_depth[dec]

            D_left = D_ends[left_edge,:]
            D_right = D_ends[right_edge,:]

            D = D_left .* D_right .* model.λ
            u0 = D

            parent_node_age = data.node_depth[anc]
            tspan = (node_age, parent_node_age)

            #prob = DifferentialEquations.ODEProblem(backward_prob, u0, tspan, pD);
            prob = DifferentialEquations.remake(prob, u0 = u0, tspan = tspan)
            sol = DifferentialEquations.solve(prob, alg, isoutofdomain = (u,p,t)->any(x->x<0,u))
            Ds[i] = sol
            sol = sol[end]
            k = sum(sol)
            sol = sol ./ k
            D_ends[i,:] = sol
            if k < 0.0
                logk = 0.0
            else
                logk = log(k)
            end
            sf[i] += logk
        end

        if verbose
            ProgressMeter.next!(prog)
        end
    end
    return(D_ends, Ds, sf, E)
end
