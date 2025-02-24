export FhBcDcModel
export eltype
#
# * Lineage-heterogeneous fossilization rate
# * Constant speciation rate
# * Constant extinction rate


struct FhBcDcModel{T1 <: Real} <: MultiStateModel
    λ::Vector{T1}
    μ::Vector{T1}
    ψ::Vector{T1}
    β::T1
end


function eltype(model::FhBcDcModel)
    return(typeof(model.β))
end


function get_speciation_rates(model::FhBcDcModel, t::Float64)
    return(model.λ)
end

function get_fossilization_rates(model::FhBcDcModel, time::Float64)
    return(model.ψ)
end

function num_parameters(model::FhBcDcModel) return 4 end

#=
function FhBcDcModel(
        λ::T1, 
        μ::T1,
        ψ::Vector{T1},
        β::T1,
    ) where {T1 <: Real}

    λv = repeat(λ, length(ψ))
    μv = repeat(μ, length(ψ))

    model = FhBcDcModel(λv, μv, ψv, β)

    return(model)
end
=#


function FhBcDc_forward_ode(du, u, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    β = model.β

    E, D = eachcol(u)
    dE, dD = eachcol(du)

    sumE = sum(E)
    sumD = sum(D)

    K = number_of_states(model)
    r = β / (K-1)

    du[:,1] .= 0.0
    du[:,2] .= 0.0

    LoopVectorization.@turbo warn_check_args=false for i in axes(u, 1)
        du[i,1] += μ[i] -(λ[i]+μ[i]+ψ[i]+β)*u[i,1] + λ[i]*u[i,1]*u[i,1] + r * (sumE - u[i,1])
        du[i,2] += +(λ[i]+μ[i]+ψ[i]+β)*u[i,2] - 2*λ[i]*u[i,2]*u[i,1] - r * (sumD - u[i,2])
    end

    nothing
end

function backward_prob(model::FhBcDcModel)
    return(FhBcDc_ode)
end


function forward_prob(model::FhBcDcModel)
    return(FhBcDc_forward_ode)
end

function extinction_prob(model::FhBcDcModel)
    return(FhBcDc_extinction_ode)
end

function FhBcDc_ode(du::Matrix{T}, u::Matrix{T}, p, t) where {T <: Real}
    model, K = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    β = model.β

    E, D = eachcol(u)
    dE, dD = eachcol(du)

    sumE = sum(E)
    sumD = sum(D)

    K = number_of_states(model)
    r = β / (K-1)

    du[:,1] .= 0.0
    du[:,2] .= 0.0

    LoopVectorization.@turbo warn_check_args=false for i in axes(u, 1)
        du[i,1] += μ[i] -(λ[i]+μ[i]+ψ[i]+β)*u[i,1] + λ[i]*u[i,1]*u[i,1] + (r * (sumE - u[i,1]))
        du[i,2] += -(λ[i]+μ[i]+ψ[i]+β)*u[i,2] + 2*λ[i]*u[i,2]*u[i,1] + (r * (sumD - u[i,2]))
    end
    
    nothing
end

function FhBcDc_extinction_ode(dE, E, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    β = model.β

    K = number_of_states(model)
    sumE = sum(E);
    r = β / (K-1);


    dE[:] = μ .- (λ.+μ.+ψ.+β).*E .+ λ.*E.*E .+ r .* (sumE .- E) 
end
