export SSEconstant, SSEdata, BDconstant

abstract type SSE <: Distributions.ContinuousUnivariateDistribution end


struct SSEconstant{T1 <: Real, T2 <: Real} <: SSE
    λ::Vector{T1}
    μ::Vector{T1}
    η::T2
end

struct SSEtimevarying <: SSE
    λ::Function
    μ::Function
    η::Function
end

struct BDconstant{T <: Real} <: Distributions.ContinuousUnivariateDistribution
    λ::T
    μ::T
end

struct SSEdata
    state_space
    trait_data
    edges::Array{Int64, 2}
    tiplab::Vector{String}
    node_depth::Vector{Float64}
    ρ::Float64
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
