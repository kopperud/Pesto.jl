

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
    α, β = x[4:5]

    μ = x[2]

    ## enforce that fossilization rate is bigger than fossilization shift rate?
    #ψmean = x[3] + β
    ψmean = x[3]

    ## enforce that 
    # * lambda hat is bigger than mu
    # * speciation rate is bigger than shift rate in speciation
    λmean = μ + x[1] + α  



    #μmean = x[1] + β*5
    #λmean = x[1] + x[2] + α*5 + β*5
    #ψmean = x[3]

    dλ = Distributions.LogNormal(log(λmean), sd)
    #dϵ = Distributions.LogNormal(log(ϵmean), sd)
    dψ = Distributions.LogNormal(log(ψmean), sd)
    
    λquantiles = Pesto.make_quantiles2(dλ, n)
    #µquantiles = Pesto.make_quantiles2(dμ, n)
    ψquantiles = Pesto.make_quantiles2(dψ, n)
    #λ, μ = allpairwise(λquantiles, µquantiles)
    #K = length(λ)
    model = FBDSconstant(λmean, μ, ψmean, λquantiles, ψquantiles, α, β)

    return(model)
end

