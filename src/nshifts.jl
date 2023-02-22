## number of shifts in the phylogeny
export compute_nshifts

function Amatrix(model, E, K, t)
    Q = -LinearAlgebra.I(K) .* model.η .+ (1 .- LinearAlgebra.I(K)) .* (model.η/(K-1))
    A = LinearAlgebra.diagm(- (model.λ .+ model.μ .- 2 .* model.λ .* E(t))) .+ Q
    return(A)
end

function Pmatrix(model, D, E, t, Δt)
    K = length(model.λ)
    A = Amatrix(model, E, K, t)

    P_unnorm = (LinearAlgebra.I(2) .- Δt .* A) .* (ones(K) * D(t)')
    #rsum = sum(I(K) .- Δt .* Amatrix(t), dims = 2)
    rsum = sum(P_unnorm, dims = 2) ## row sum
    P = P_unnorm ./ rsum
    return(P)
end

function compute_nshifts(model, data, Ds, Ss; ntimes = 200)
    E = extinction_probability(model, data)
    nbranches = size(data.edges)[1]
    K = length(model.λ)
    nshifts = zeros(nbranches)

    ProgressMeter.@showprogress for edge_idx in 1:nbranches
        a = Ds[edge_idx].t[end]
        b = Ds[edge_idx].t[1]

        times = collect(range(a, b, length = ntimes))
        Δt = times[2] - times[1]

        nshift = 0.0
        Ps = zeros(ntimes, K, K)
        for i in 1:(ntimes-1)
            #P = trans_prob8(model, Ds[edge_idx], times[i], Δt)
            P = Pmatrix(model, Ds[edge_idx], E, times[i], Δt)
            state_prob = Ss[edge_idx](times[i])

            ## this one is good
            #nshift += abs((P' * state_prob)[1] - (P' * state_prob)[2])

            P0 = (1 .- LinearAlgebra.I(K)) .* P
            L = LinearAlgebra.LowerTriangular(P0)
            U = LinearAlgebra.UpperTriangular(P0)

            L1 = L .* (state_prob * ones(K)')
            U1 = U .* (state_prob * ones(K)')

            nshift += sum(abs.(L1 .- U1'))
        end
        nshifts[edge_idx] = nshift
    end

    ## reorder to ape node indices
    ancestors = Diversification.make_ancestors(data)

    node_nshifts = zeros(maximum(data.edges))
    for i in 1:maximum(data.edges)
        if i == length(data.tiplab)+1
            node_nshifts[i] = 0.0
        else
            edge_idx = ancestors[i]
            node_val = nshifts[edge_idx]
            node_nshifts[i] = node_val
        end
    end

    return(node_nshifts)
end