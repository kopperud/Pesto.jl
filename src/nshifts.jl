## number of shifts in the phylogeny
export compute_nshifts
export state_shifts

function Amatrix(model, E, K, t)
    Q = -LinearAlgebra.I(K) .* model.η .+ (1 .- LinearAlgebra.I(K)) .* (model.η/(K-1))
    A = LinearAlgebra.diagm(- (model.λ .+ model.μ) .+ 2 .* model.λ .* E(t)) .+ Q
    return(A)
end

function Pmatrix(model, D, E, t, Δt)
    K = length(model.λ)
    A = Amatrix(model, E, K, t)

    P_unnorm = (LinearAlgebra.I(K) .- Δt .* A) .* (ones(K) * D(t)')
    csum = sum(P_unnorm, dims = 1)' * ones(K)' ## column sum
    P = P_unnorm ./ csum
    return(P)
end

function state_shifts(model, data; ntimeslices = 500, ape_order = true)
    Ds, Fs = backwards_forwards_pass(model, data);
    Ss = ancestral_state_probabilities(data, Ds, Fs);

    state_shifts(model, data, Ds, Ss; ntimeslices = ntimeslices, ape_order = ape_order)
end

function state_shifts(model, data, Ds, Ss; ntimeslices = 500, ape_order = true)
    E = extinction_probability(model, data)
    nbranches = size(data.edges)[1]
    K = length(model.λ)
    nshifts = zeros(nbranches, K, K)
    th = maximum(data.branching_times)

    for edge_idx in 1:nbranches
        a = Ds[edge_idx].t[end]
        b = Ds[edge_idx].t[1]

        branch_ntimeslices = Int64(round(ntimeslices * data.branch_lengths[edge_idx] / th, RoundUp))

        times = collect(range(a, b, length = branch_ntimeslices+2))
        ntimes = length(times)
        Δt = times[2] - times[1]

        nshift = zeros(K, K)
        for i in 1:(ntimes-1)
            P = Pmatrix(model, Ds[edge_idx], E, times[i], Δt)
            state_prob = Ss[edge_idx](times[i]+Δt)
            
            W = state_prob * ones(K)'
            Ntilde = (1.0 .- LinearAlgebra.I(K)) .* P .* W
            Nhat = transpose(Ntilde)
            
            nshift[:,:] += Nhat
        end
        nshifts[edge_idx,:,:] = nshift
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

function compute_nshifts(model, data; ntimeslices = 500, ape_order = true)
    Ds, Fs = backwards_forwards_pass(model, data);
    Ss = ancestral_state_probabilities(data, Ds, Fs);

    compute_nshifts(model, data, Ds, Ss; ntimeslices = ntimeslices, ape_order = ape_order)
end

function compute_nshifts(model, data, Ds, Ss; ntimeslices = 500, ape_order = true)
    nshifts = state_shifts(model, data, Ds, Ss; ntimeslices = ntimeslices, ape_order = ape_order)
    res = sum(nshifts, dims = 2:3)[:,1,1]
    return(res)
end