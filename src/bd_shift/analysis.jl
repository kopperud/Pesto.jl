## wrapper for analysis
export pesto
export pesto_twostep

export empirical_bayes

function empirical_bayes(data::SSEdata; n = 6, sd = 0.587)
    λml, μml = estimate_constant_bdp(data)

    dλ = Distributions.LogNormal(log(λml), sd)
    dμ = Distributions.LogNormal(log(μml), sd)

    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)
    η = optimize_eta(λ, µ, data)
    model = BDSconstant(λ, μ, η)
    return(model)
end

function pesto_twostep(data::SSEdata; n = 6, sd = 0.587)
    model = empirical_bayes(data; n = n)

    rates = birth_death_shift(model, data);
    return(model, rates)
end

function pesto(data::SSEdata; n = 6, sd = 0.587)
    optres, model, i = optimize_hyperparameters(data; n = n, sd = sd, n_attempts = 5)

    rates = birth_death_shift(model, data);
    return(model, rates)
end
