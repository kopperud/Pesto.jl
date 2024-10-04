export Node
export Tip 

export Branch
export Root
export InternalNode

abstract type AbstractNode end
abstract type AbstractBranch end
abstract type InternalNode <: AbstractNode end

mutable struct Branch <: AbstractBranch
    index::Int64
    inbounds::AbstractNode
    outbounds::AbstractNode

    time::Float64
    Branch() = new()
end

mutable struct Tip <: AbstractNode
    index::Int64
    inbounds::Branch
    label::String
    sampling_probability::Float64
    is_fossil::Bool

    Tip() = new()
end


mutable struct Node <: InternalNode
    index::Int64
    inbounds::Branch
    children::Vector{Branch}
    Node() = new()
end

mutable struct Root <: InternalNode
    index::Int64
    children::Vector{Branch}
    Root() = new()
end
