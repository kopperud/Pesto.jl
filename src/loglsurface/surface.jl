export logl_slice

function lognormal_pairwise(λmean, μmean, η; n = 6, sd = 0.587)
    dλ = Distributions.LogNormal(log(λmean), sd)
    dμ = Distributions.LogNormal(log(μmean), sd)
    λq = make_quantiles(dλ, n)
    μq = make_quantiles(dμ, n)

    λ, μ = allpairwise(λq, μq)

    model = SSEconstant(λ, μ, η)

    return(model)
end

function logl_slice(
        λmin::Float64, 
        λmax::Float64,
        μmin::Float64, 
        μmax::Float64,
        η::Float64,
        data::SSEdata;
        i = 10,
        j = 15,
        n = 6,
        sd = 0.587
    )

    λs = collect(range(λmin, λmax; length = i))
    μs = collect(range(μmin, μmax; length = i))

    #prog = ProgressMeter.Progress(n*m)A

    logls = zeros(i, j)

    for ii in 1:i
        for jj in 1:j
             
            λ = λs[ii]
            μ = μs[ii]

            model = lognormal_pairwise(λ, μ, η; n = n, sd = sd)
            Threads.@spawn logls[ii,jj] = logL_root(model, data, multithread = false)
    #        Progress.next!(prog)
        end
    end

    #Progress.Meterfinish!(prog)

    return(logls)

end
