export SSEconstant, SSEdata, BDconstant
export BDSconstant, FBDSconstant

#abstract type SSE <: Distributions.ContinuousUnivariateDistribution end
abstract type Model end

abstract type ConstantModel <: Model end
abstract type TimevaryingModel <: Model end

#struct BDconstant{T <: Real}#<: Distributions.ContinuousUnivariateDistribution
struct BDconstant{T <: Real} <: ConstantModel
    λ::T
    μ::T
end

struct BDSconstant{T1 <: Real, T2 <: Real} <: ConstantModel
    λ::Vector{T1}
    μ::Vector{T1}
    η::T2
end

struct BDStimevarying <: TimevaryingModel
    λ::Function
    μ::Function
    η::Function
end
#=
struct FBDSconstant{T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real} <: ConstantModel
    λ::Vector{T1}
    μ::Vector{T2}
    ψ::Vector{T3}
    η::T4
end
=#

struct FBDSconstant{T1 <: Real, T2 <: Real, T3 <: Real} <: ConstantModel
    λmean::T1
    ψmean::T1
    #μmean::T1
    λ::Vector{T2}
    μ::Vector{T2}
    ψ::Vector{T2}
    α::T3
    β::T3
    #γ::T3
    Q::SparseArrays.SparseMatrixCSC{T3, Int64}
    Qα::SparseArrays.SparseMatrixCSC{Int64, Int64}
    Qβ::SparseArrays.SparseMatrixCSC{Int64, Int64}
    #Qγ::SparseArrays.SparseMatrixCSC{Int64, Int64}
end

function FBDSconstant(
        λmean::T1,
        μ::T1,
        ψmean::T1,
        λ::Vector{T2},
        ψ::Vector{T2},
        α::T3,
        β::T3,
        #γ::T3,
    ) where {T1 <: Real, T2 <: Real, T3 <: Real}
    
    #λ = r ./ (1 - ϵ)
    #μ = λ .- r

    λv, ψv = allpairwise(λ, ψ)
    μv = repeat([μ], length(λv))
    #μv, _ = allpairwise(μ, ψ)

    Q, Qα, Qβ = Qmatrix(λ, ψ, α, β)

    model = FBDSconstant(λmean, ψmean, λv, μv, ψv, α, β, Q, Qα, Qβ)

    return(model)
end

## type aliases to work with the old names
const SSE = Model
const SSEconstant = BDSconstant
const SSEtimevarying = BDStimevarying


struct SSEdata
    state_space
    trait_data
    edges::Array{Int64, 2}
    tiplab::Vector{String}
    node_depth::Vector{Float64}
    sampling_probability::Float64
    branch_lengths::Vector{Float64}
    branching_times::Vector{Float64}
    po::Vector{Int64}
    Nnode::Int64
end

struct phylo
    edge::Matrix{Int64}
    edge_length::Vector{Float64}
    Nnode::Int64
    tip_label::Vector{String}
    node_depths::Vector{Float64}
    branching_times::Vector{Float64}
    po::Vector{Int64}
end

struct SSEresult
    phy
    rates
end
