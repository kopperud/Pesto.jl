

function newmodel(x::Vector{T}; n = 6, sd = 0.587) where {T <: Real}
    η = x[1]
    μmean = x[2]
    λmean = maximum([5*x[1], x[2]]) + x[3]
            
    dλ = Distributions.LogNormal(log(λmean), sd)
    dμ = Distributions.LogNormal(log(μmean), sd)
    
    λquantiles = Pesto.make_quantiles2(dλ, n)
    µquantiles = Pesto.make_quantiles2(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)
    model = BDSconstant(λ, μ, η)

    return(model)
end

function newmodel_fbd(x::Vector{T}; n = 6, sd = 0.587) where {T <: Real}
    η = x[1]
    μmean = x[2]
    λmean = maximum([5*x[1], x[2]]) + x[3]
    ψ = x[4]
    
    dλ = Distributions.LogNormal(log(λmean), sd)
    dμ = Distributions.LogNormal(log(μmean), sd)
    
    λquantiles = Pesto.make_quantiles2(dλ, n)
    µquantiles = Pesto.make_quantiles2(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)
    K = length(λ)
    model = FBDSconstant(λ, μ, repeat([ψ], K), η)

    return(model)
end

