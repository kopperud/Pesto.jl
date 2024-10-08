

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
    α, β, γ = x[4:6]
    μmean = x[1] + β*5
    #λmean = maximum([α+β, x[1]]) + x[2]
    λmean = x[1] + x[2] + α*5 + β*5
    ψmean = x[3]

    #=
    η = x[1]
    μmean = x[2]
    λmean = maximum([5*x[1], x[2]]) + x[3]
    ψ = x[4]
    =#
    
    dλ = Distributions.LogNormal(log(λmean), sd)
    dμ = Distributions.LogNormal(log(μmean), sd)
    dψ = Distributions.LogNormal(log(ψmean), sd)
    
    λquantiles = Pesto.make_quantiles2(dλ, n)
    µquantiles = Pesto.make_quantiles2(dμ, n)
    ψquantiles = Pesto.make_quantiles2(dψ, n)
    #λ, μ = allpairwise(λquantiles, µquantiles)
    #K = length(λ)
    model = FBDSconstant(λmean, μmean, ψmean, λquantiles, μquantiles, ψquantiles, α, β, γ)

    return(model)
end

