## wrapper foranalysis
export pesto

function pesto(data; n = 6, sd = 0.587)
    λml, μml = estimate_constant_bdp(data)

    dλ = Distributions.LogNormal(log(λml), sd)
    dμ = Distributions.LogNormal(log(μml), sd)

    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)
    η = optimize_eta(λ, µ, data)
    model = SSEconstant(λ, μ, η)

    rates = birth_death_shift(model, data);
    return(model, rates)
end