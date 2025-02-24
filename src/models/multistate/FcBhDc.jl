export FcBhDcModel
export eltype
#
# * Constant fossilization rate
# * Lineage-heterogeneous speciation rate
# * Constant extinction rate


struct FcBhDcModel{T1 <: Real} <: MultiStateModel
    λ::Vector{T1}
    μ::Vector{T1}
    ψ::Vector{T1}
    α::T1
end


function eltype(model::FcBhDcModel)
    return(typeof(model.α))
end


function get_speciation_rates(model::FcBhDcModel, t::Float64)
    return(model.λ)
end

function get_fossilization_rates(model::FcBhDcModel, time::Float64)
    return(model.ψ)
end

function num_parameters(model::FcBhDcModel) return 4 end

function FcBhDc_forward_ode(du, u, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    α = model.α

    E, D = eachcol(u)
    dE, dD = eachcol(du)

    sumE = sum(E)
    sumD = sum(D)

    K = number_of_states(model)
    r = α / (K-1)

    du[:,1] .= 0.0
    du[:,2] .= 0.0

    LoopVectorization.@turbo warn_check_args=false for i in axes(u, 1)
        du[i,1] += μ[i] -(λ[i]+μ[i]+ψ[i]+α)*u[i,1] + λ[i]*u[i,1]*u[i,1] + r * (sumE - u[i,1])
        du[i,2] += +(λ[i]+μ[i]+ψ[i]+α)*u[i,2] - 2*λ[i]*u[i,2]*u[i,1] - r * (sumD - u[i,2])
    end

    nothing
end

function backward_prob(model::FcBhDcModel)
    return(FcBhDc_ode)
end


function forward_prob(model::FcBhDcModel)
    return(FcBhDc_forward_ode)
end

function extinction_prob(model::FcBhDcModel)
    return(FcBhDc_extinction_ode)
end

function FcBhDc_ode(du::Matrix{T}, u::Matrix{T}, p, t) where {T <: Real}
    model, K = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    α = model.α

    E, D = eachcol(u)
    dE, dD = eachcol(du)

    sumE = sum(E)
    sumD = sum(D)

    K = number_of_states(model)
    r = α / (K-1)

    du[:,1] .= 0.0
    du[:,2] .= 0.0

    LoopVectorization.@turbo warn_check_args=false for i in axes(u, 1)
        du[i,1] += μ[i] -(λ[i]+μ[i]+ψ[i]+α)*u[i,1] + λ[i]*u[i,1]*u[i,1] + (r * (sumE - u[i,1]))
        du[i,2] += -(λ[i]+μ[i]+ψ[i]+α)*u[i,2] + 2*λ[i]*u[i,2]*u[i,1] + (r * (sumD - u[i,2]))
    end
    
    nothing
end

function FcBhDc_extinction_ode(dE, E, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    α = model.α

    K = number_of_states(model)
    sumE = sum(E);
    r = α / (K-1);


    dE[:] = μ .- (λ.+μ.+ψ.+α).*E .+ λ.*E.*E .+ r .* (sumE .- E) 
end

