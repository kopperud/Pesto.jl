export makebins

function makebins(N0, model, lower, upper; filter = "", nbins = 18)
    λ = model.λ
    μ = model.μ
    K = length(λ)
    N = deepcopy(N0)

    ## sp
    Δλ = λ * ones(K)' .- ones(K) * λ'
    ## ex
    Δμ = μ * ones(K)' .- ones(K) * μ'
    ## netdiv
    r = λ .- μ
    Δr = r * ones(K)' .- ones(K) * r'
    ## relext
    ϵ = μ ./ λ
    Δϵ = ϵ * ones(K)' .- ones(K) * ϵ'

    ## remove diagonals just to be sure
    for i in 1:K
        N[i,i] = 0.0
    end

    if filter == "speciation"
        is_zero = Δλ .!= 0
        N[is_zero] .= 0
    elseif filter == "extinction"
        is_zero = Δμ .!= 0
        N[is_zero] .= 0
    elseif filter == "speciation+extinction"
        is_zero_ex = Δμ .!= 0
        is_zero_sp = Δλ .!= 0

        is_zero = is_zero_ex .& is_zero_sp
        N[is_zero] .= 0
    end

    borders = collect(range(lower, upper; length = nbins+1))
    mids = [(borders[i]+borders[i+1])/2 for i in 1:nbins]
    bins = zeros(nbins, 4)

    for (j, Δx) in enumerate([Δλ, Δμ, Δr, Δϵ])
        for i in 1:nbins
            in_bin = (Δx .> borders[i]) .& (Δx .<= borders[i+1])
            not_zero = N .> 0
            idx = in_bin .& not_zero
            bins[i,j] = sum(N[idx])
        end
    end
    return (mids, bins)
end
