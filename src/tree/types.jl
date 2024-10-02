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
    species_name::String
    is_fossil::Bool

    Tip() = new()
end

mutable struct Node <: InternalNode
    index::Int64
    inbounds::Branch
    left::Branch
    right::Branch
    Node() = new()
end

mutable struct Root <: InternalNode
    index::Int64
    left::Branch
    right::Branch
    Root() = new()
end
