export FcBhDhModel
#
# * Constant fossilization rate
# * Lineage-heterogeneous speciation rate
# * Lineage-heterogeneous extinction rate


struct FcBhDhModel{T1<:Real,T2<:Real,T3<:Real} <: MultiStateModel
    λ::Vector{T1}
    μ::Vector{T1}
    ψ::Vector{T1}
    α::T2
    β::T2
    Q::Matrix{T3}
end

const FBDSconstant = FcBhDhModel

function eltype(model::FcBhDhModel)
    return (typeof(model.α))
end

function get_speciation_rates(model::FcBhDhModel, t::Float64)
    return (model.λ)
end

function get_fossilization_rates(model::FcBhDhModel, time::Float64)
    return (model.ψ)
end

function num_parameters(model::FcBhDhModel)
    return 5
end

function FcBhDhModel(
    λ::Vector{T1}, ## these are of length n, not n*n or n^3
    μ::Vector{T1},
    ψ::Vector{T1},
    α::T2,
    β::T2,
) where {T1<:Real,T2<:Real}

    #λ = r ./ (1 - ϵ)
    #μ = λ .- r

    λv, μv = allpairwise(λ, μ)
    μv = repeat(μ, length(μ))
    ψv = repeat(ψ, length(ψ))
    #μv, _ = allpairwise(μ, ψ)

    n = length(λ)
    Q = Qmatrix(α, β, n)

    model = FcBhDhModel(λv, μv, ψv, α, β, Q)

    return (model)
end


function backward_prob(model::FcBhDhModel)
    return (FhBhDc_ode)
end

function forward_fossil_ode(du, u, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    Q = model.Q

    E, F = eachcol(u)
    dE, dF = eachcol(du)

    fastmv!(dE, Q, .- E)
    fastmv!(dF, Q, .- F)

    LoopVectorization.@turbo warn_check_args=false for i in axes(u, 1)
        du[i, 1] += - μ[i] + (λ[i]+μ[i]+ψ[i])*u[i, 1] - λ[i]*u[i, 1]*u[i, 1]
        du[i, 2] += +(λ[i]+μ[i]+ψ[i])*u[i, 2] - 2*λ[i]*u[i, 2]*u[i, 1]
    end

    nothing
end


function forward_prob(model::FcBhDhModel)
    return (forward_fossil_ode)
end

function extinction_prob(model::FcBhDhModel)
    return (extinction_fossil_ode)
end

function FhBhDc_ode(du::Matrix{T}, u::Matrix{T}, p, t) where {T<:Real}
    model, K = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    #η = model.η
    Q = model.Q

    #=
    Et = E(t)
    #dD[:] .= - (λ .+ μ .+ ψ .+ η) .* D .+ 2 .* λ .* D .* Et .+ (η/(K-1)) .* (sum(D) .- D)
    dD[:] .= - (λ .+ μ .+ ψ) .* D .+ 2 .* λ .* D .* Et .+ Q * D
    =#

    E, D = eachcol(u)
    dE, dD = eachcol(du)

    fastmv!(dE, Q, E)
    fastmv!(dD, Q, D)

    LoopVectorization.@turbo warn_check_args=false for i in axes(u, 1)
        du[i, 1] += μ[i] - (λ[i]+μ[i]+ψ[i])*u[i, 1] + λ[i]*u[i, 1]*u[i, 1]
        du[i, 2] += -(λ[i]+μ[i]+ψ[i])*u[i, 2] + 2*λ[i]*u[i, 2]*u[i, 1]
    end

    nothing
end

#=
function backward_fossil2_ode(dD, D, p, t)
    model, K, E = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    Q = model.Q

    Et = E(t)
    dD[:] = - (λ .+ μ .+ ψ) .* D .+ 2 .* λ .* D .* Et .+ Q * D 
end
=#


function extinction_fossil_ode(dE, E, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    Q = model.Q
    K = number_of_states(model)

    dE[:] = μ .- (λ .+ μ .+ ψ) .* E .+ λ .* E .* E .+ Q * E
end

