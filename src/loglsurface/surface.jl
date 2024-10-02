export logl_slice

function lognormal_pairwise(λmean, μmean, η; n = 6, sd = 0.587)
    dλ = Distributions.LogNormal(log(λmean), sd)
    dμ = Distributions.LogNormal(log(μmean), sd)
    λq = make_quantiles(dλ, n)
    μq = make_quantiles(dμ, n)

    λ, μ = allpairwise(λq, μq)

    model = BDSconstant(λ, μ, η)

    return(model)
end

export plot_loglsurface

function plot_loglsurface() end

function logl_slice(
        λ::Vector{Float64},
        μ::Vector{Float64},
        η::Vector{Float64},
        data::SSEdata;
        n_λ = 10,
        n_μ = 11,
        n_η = 12,
        n_categories = 6,
        sd = 0.587
    )
    @assert length(λ) == 3
    @assert length(μ) == 3
    @assert length(η) == 3

    λ_fixed, λ_min, λ_max = λ
    μ_fixed, μ_min, μ_max = μ
    η_fixed, η_min, η_max = η

    @assert λ_min < λ_fixed < λ_max
    @assert μ_min < μ_fixed < μ_max
    @assert η_min < η_fixed < η_max

    λs = collect(range(λ_min, λ_max; length = n_λ))
    μs = collect(range(μ_min, μ_max; length = n_μ))
    ηs = collect(range(η_min, η_max; length = n_η))

    prog = ProgressMeter.Progress(n_λ * n_μ * n_η; desc =  "hello")

    logls_fixed_lambda = zeros(n_μ, n_η)
    logls_fixed_mu = zeros(n_λ, n_η)
    logls_fixed_eta = zeros(n_λ, n_μ)

    ## fixed lambda
    Threads.@sync for i in 1:n_μ, j in 1:n_η
        Threads.@spawn begin
            μi = μs[i]
            ηj = ηs[j]

            model = lognormal_pairwise(λ_fixed, μi, ηj; n = n_categories, sd = sd)
            logls_fixed_lambda[i, j] = logL_root(model, data, multithread = false)
            ProgressMeter.next!(prog)
        end
    end

    ## fixed μ
    Threads.@sync for i in 1:n_λ, j in 1:n_η
        Threads.@spawn begin
            λi = λs[i]
            ηj = ηs[j]

            model = lognormal_pairwise(λi, μ_fixed, ηj; n = n_categories, sd = sd)
            logls_fixed_mu[i, j] = logL_root(model, data, multithread = false)
            ProgressMeter.next!(prog)
        end
    end

    ## fixed η
    Threads.@sync for i in 1:n_λ, j in 1:n_μ
        Threads.@spawn begin
            λi = λs[i]
            μj = μs[j]

            model = lognormal_pairwise(λi, μj, η_fixed; n = n_categories, sd = sd)
            logls_fixed_eta[i,j] = logL_root(model, data, multithread = false)
            ProgressMeter.next!(prog)
        end
    end

    ProgressMeter.finish!(prog)

    logls = [
             logls_fixed_lambda,
             logls_fixed_mu,
             logls_fixed_eta
            ]

    vectors = [
               λs, μs, ηs
              ]

    return((vectors, logls))

end
