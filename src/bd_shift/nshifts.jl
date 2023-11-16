## number of shifts in the phylogeny
export compute_nshifts
export state_shifts

function state_shifts(model::SSE, data::SSEdata; ape_order = true)
    Ds, Fs = backwards_forwards_pass(model, data);
    Ss = ancestral_state_probabilities(data, Ds, Fs);

    state_shifts(model, data, Ds, Ss; ape_order = ape_order)
end

function state_shifts(model::SSE, data::SSEdata, Ds, Ss; alg = OrdinaryDiffEq.Tsit5(), ape_order = true)
    nbranches = size(data.edges)[1]
#    K = length(model.λ)
    K = number_of_states(model)    
    nshifts = zeros(nbranches, K, K)
    ode = shift_problem(model)

    #Threads.@threads for edge_idx in 1:nbranches
    for edge_idx in 1:nbranches
        a = Ds[edge_idx].t[end]
        b = Ds[edge_idx].t[1]
        tspan = (a,b)

        N0 = zeros(K,K)
        p = (model.η, K, Ss[edge_idx], Ds[edge_idx])

        prob = OrdinaryDiffEq.ODEProblem(ode, N0, tspan, p)
        sol = OrdinaryDiffEq.solve(prob, alg, isoutofdomain = (u,p,t)->any(x->x<0,u))

        nshifts[edge_idx,:,:] = sol[end]
    end

    if ape_order
        ## reorder to ape node indices
        ancestors = make_ancestors(data)

        node_nshifts = zeros(maximum(data.edges), K, K)
        for i in 1:maximum(data.edges)
            if i == length(data.tiplab)+1
                node_nshifts[i,:,:] .= 0.0
            else
                edge_idx = ancestors[i]
                node_val = nshifts[edge_idx,:,:]
                node_nshifts[i,:,:] = node_val
            end
        end
        return(node_nshifts)
    else
        return(nshifts)
    end
end

function compute_nshifts(model, data; ape_order = true)
    Ds, Fs = backwards_forwards_pass(model, data);
    Ss = ancestral_state_probabilities(data, Ds, Fs);

    compute_nshifts(model, data, Ds, Ss; ape_order = ape_order)
end

function compute_nshifts(model, data, Ds, Ss; ape_order = true)
    nshifts = state_shifts(model, data, Ds, Ss; ape_order = ape_order)
    res = sum(nshifts, dims = 2:3)[:,1,1]
    return(res)
end
