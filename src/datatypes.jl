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

struct FBDSconstant{T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real} <: ConstantModel
    λ::Vector{T1}
    μ::Vector{T2}
    ψ::Vector{T3}
    η::T4
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
