export SSEdata 

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

#=
struct BDSconstantQ{T1 <: Real, T2 <: Real} <: ConstantModel
    λ::Vector{T1}
    μ::Vector{T1}
    Q::Matrix{T2}
end
=#

#=
struct BDStimevarying <: TimevaryingModel
    λ::Function
    μ::Function
    η::Function
end
=#

#=
struct FBDSconstant{T1 <: Real, T2 <: Real, T3 <: Real, T4 <: Real} <: ConstantModel
    λ::Vector{T1}
    μ::Vector{T2}
    ψ::Vector{T3}
    η::T4
end
=#

## type aliases to work with the old names
#const SSEtimevarying = BDStimevarying



