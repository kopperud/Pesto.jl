export SSEconstant, SSEdata, BDconstant

abstract type SSE <: Distributions.ContinuousUnivariateDistribution end


struct SSEconstant <: SSE
    λ
    μ
    η
end

struct SSEtimevarying <: SSE
    λ::Function
    μ::Function
    η::Function
end

struct BDconstant <: Distributions.ContinuousUnivariateDistribution
    λ
    μ
end

struct SSEdata
    state_space
    trait_data
    edges
    tiplab
    node_depth
    ρ
    branch_lengths
    branching_times
    po
    Nnode
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
