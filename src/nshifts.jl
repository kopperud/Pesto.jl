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

    P_unnorm = (LinearAlgebra.I(K) .- Δt .* A) .* (D(t) * ones(K)')
    #rsum = sum(I(K) .- Δt .* Amatrix(t), dims = 2)
    rsum = sum(P_unnorm, dims = 1)' * ones(K)' ## row sum
    #rsum = sum(P_unnorm, dims = 1)
    P = P_unnorm' ./ rsum
    return(P)
end

function compute_nshifts(model, data, Ds, Ss; ntimes = 100, ape_order = true)
    E = extinction_probability(model, data)
    nbranches = size(data.edges)[1]
    K = length(model.λ)
    nshifts = zeros(nbranches)

    for edge_idx in 1:nbranches
        a = Ds[edge_idx].t[end]
        b = Ds[edge_idx].t[1]

        times = collect(range(a, b, length = ntimes))
        Δt = times[2] - times[1]

        nshift = 0.0
        for i in 1:(ntimes-1)
            P = Pmatrix(model, Ds[edge_idx], E, times[i], Δt)
            state_prob = Ss[edge_idx](times[i])

            P0 = (1 .- LinearAlgebra.I(K)) .* P
            ## this one is good
            #nshift += abs((P0' * state_prob)[1] - (P0' * state_prob)[2])

            W = state_prob * ones(K)'
            #W = ones(K) * state_prob'

            P1 = P0 .* W

            L = LinearAlgebra.LowerTriangular(P1)
            U = LinearAlgebra.UpperTriangular(P1)
            
            #W = ones(K) * state_prob'
            #W = ones(K) * ones(K)'
            
            #L1 = L .* W
            #U1 = U .* W

            nshift += sum(abs.(L .+ U'))
            #nshift += sum(abs.(U1))
        end
        nshifts[edge_idx] = nshift
    end
    println("asd")

    if ape_order
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
    else
        return(nshifts)
    end
end