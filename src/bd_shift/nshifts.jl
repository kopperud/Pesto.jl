## number of shifts in the phylogeny
export compute_nshifts
export state_shifts
export state_shifts_simple

function state_shifts(model::SSE, data::SSEdata)
    Ds, Fs = backwards_forwards_pass(model, data);
    #Ss = ancestral_state_probabilities(data, Ds, Fs);

    state_shifts(model, data, Ds, Fs)
end

function state_shifts_simple(model::SSE, data::SSEdata, Ds, Fs; alg = OrdinaryDiffEq.Tsit5())
    nbranches = size(data.edges)[1]
    K = number_of_states(model)    
    nshifts = zeros(nbranches)
    ode = shift_problem_simple(model)

    Threads.@threads for edge_idx in 1:nbranches
    #for edge_idx in 1:nbranches
        a = Ds[edge_idx].t[end]
        b = Ds[edge_idx].t[1]
        tspan = (a,b)

        N0 = [0.0]

        p = (model.η, K, Ds[edge_idx], Fs[edge_idx]);

        prob = OrdinaryDiffEq.ODEProblem(ode, N0, tspan, p);
        sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg)

        nshifts[edge_idx] = sol[end][1]
    end

    return(nshifts)
end


function reorder_ape(nshifts::Array{Float64, 1}, data::SSEdata)
    ancestors = make_ancestors(data)

    n_nodes = data.Nnode + length(data.tiplab)
    root_index = length(data.tiplab)+1
    #n_nodes = maximum(data.edges)

    node_nshifts = zeros(Float64, n_nodes)
    for i in 1:n_nodes
        if i == root_index
            node_nshifts[i] = 0.0
        else
            edge_idx = ancestors[i]
            node_val = nshifts[edge_idx]
            node_nshifts[i] = node_val
        end
    end
    return(node_nshifts)
end

function reorder_ape(nshifts::Array{Float64, 3}, data::SSEdata, K::Int64)
    ancestors = make_ancestors(data)

    n_nodes = data.Nnode + length(data.tiplab)
    root_index = length(data.tiplab)+1
    #n_nodes = maximum(data.edges)

    node_nshifts = zeros(Float64, n_nodes, K, K)
    for i in 1:n_nodes
        if i == root_index
            node_nshifts[i,:,:] .= 0.0
        else
            edge_idx = ancestors[i]
            node_val = nshifts[edge_idx,:,:]
            node_nshifts[i,:,:] = node_val
        end
    end
    return(node_nshifts)
end

function state_shifts(model::SSE, data::SSEdata, Ds, Fs; alg = OrdinaryDiffEq.Tsit5())
    nbranches = size(data.edges)[1]
    K = number_of_states(model)    
    nshifts = zeros(Float64, nbranches, K, K)
    ode = shift_problem(model)

    Threads.@threads for edge_idx in 1:nbranches
        a = Ds[edge_idx].t[end]
        b = Ds[edge_idx].t[1]
        tspan = (a,b)
        N0 = zeros(K,K)

        p = (model.η, K, Ds[edge_idx], Fs[edge_idx])

        prob = OrdinaryDiffEq.ODEProblem(ode, N0, tspan, p)
        sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg)

        nshifts[edge_idx,:,:] = sol[end]
    end

    return(nshifts)
end

function compute_nshifts(model::SSE, data::SSEdata, Ds, Fs)
    nshifts = state_shifts_simple(model, data, Ds, Fs)
    return(nshifts)
end


function nshifts_simple_whole_branch(
        model::SSE, 
        data::SSEdata, 
        Ds, 
        Fs; 
        alg = OrdinaryDiffEq.Tsit5())

    nbranches = size(data.edges)[1]
    K = number_of_states(model)    
    nshifts = Dict() ## problem: not type stable
    ode = shift_problem_simple(model)

    for edge_idx in 1:nbranches
        a = Ds[edge_idx].t[end]
        b = Ds[edge_idx].t[1]
        tspan = (a,b)

        N0 = Float64[0.0] 

        p = (model.η, K, Ds[edge_idx], Fs[edge_idx]);

        prob = OrdinaryDiffEq.ODEProblem(ode, N0, tspan, p);
        sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = notneg)

        nshifts[edge_idx] = sol
    end

    return(nshifts)
end
