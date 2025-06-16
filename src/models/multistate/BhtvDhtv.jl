export BhtvDhtvModel, BDSconstant, SSEconstant

struct BhtvDhtvModel <: MultiStateModel
    λ::Function
    μ::Function
    η::Function
end

function eltype(model::BhtvDhtvModel)
    return(typeof(model.λ(0.0)[1]))
end

function get_speciation_rates(model::BhtvDhtvModel, t::Float64)
    return(model.λ(t))
end

function get_extinction_rates(model::BhtvDhtvModel, t::Float64)
    return(model.μ(t))
end

function get_fossilization_rates(model::BhtvDhtvModel, time::Float64)
    error("can not get fossilization rate for a birth-death-shift model. either use an FBD model or don't include fossils in the tree.")
    #return(model.ψ)
end

function BhtvDhtv_extinction_ode(dE, E, p, t)
    model, K = p
    λ = model.λ(t)
    μ = model.μ(t)
    η = model.η(t)

    sumE = sum(E)

    LoopVectorization.@turbo warn_check_args=false for i in 1:K
        dE[i] = μ[i] - (λ[i] + μ[i] + η) * E[i] + λ[i] * E[i] * E[i] + (η/(K-1)) * (sumE - E[i]) 
    end
    nothing
end

function extinction_prob(model::BhtvDhtvModel)
    return(BhtvDhtv_extinction_ode)
end


## Probability of of observing the branch at time `t`
## * We solve this equation in the postorder traversal
function BhtvDhtv_backward_ode(
        du::Matrix{T}, 
        u::Matrix{T}, 
        p, 
        t::Float64
    ) where {T <: Real}

    model, K = p
    E, D = eachcol(u)
    sumE = sum(E)
    sumD = sum(D)

    λ = model.λ(t)
    μ = model.μ(t)
    η = model.η(t)

    r = η / (K-1)

    LoopVectorization.@turbo warn_check_args=false for i in axes(du, 1)
        du[i,1] = μ[i] -(λ[i]+μ[i]+η)*u[i,1] + λ[i]*u[i,1]*u[i,1] + r*(sumE-u[i,1]) 
        du[i,2] = -(λ[i]+μ[i]+η)*u[i,2] + 2*λ[i]*u[i,2]*u[i,1] + r*(sumD-u[i,2])
    end

    nothing
end

function backward_prob(model::BhtvDhtvModel)
    return(BhtvDhtv_backward_ode)
end

## This ODE is the previous one times minus one
## * We solve this equation in the preorder traversal, albeit with different starting values for each branch
function BhtvDhtv_forward_ode(
        dF::Vector{T}, 
        F::Vector{T}, 
        p, 
        t::Float64) where {T <: Real}

    model, E, K = p
    λ = model.λ(t)
    μ = model.μ(t)
    η = model.η(t)

    Et = E(t)

    sumF = sum(F)

    r = η / (K-1)

    LoopVectorization.@turbo warn_check_args=false for i in eachindex(dF)
        #du[i,1] = μ[i] - (λ[i] + μ[i] + η) * u[i,1] + λ[i] * u[i,1] * u[i,1] + r * (sumE - u[i,1]) 
        ## F
        dF[i] = +(λ[i]+μ[i]+η)*F[i] - 2*λ[i]*F[i]*Et[i] - r *(sumF - F[i])
    end
    nothing
end

function forward_prob(model::BhtvDhtvModel)
    return(BhtvDhtv_forward_ode)
end

