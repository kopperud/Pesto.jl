export Node
export Tip 

export Branch
export Root
export InternalNode

abstract type AbstractNode end
abstract type AbstractBranch end

########################################################
##
##              branch
##
########################################################
mutable struct Branch <: AbstractBranch
    index::Int64
    inbounds::AbstractNode
    outbounds::AbstractNode

    time::Float64
    Branch() = new()
end

########################################################
##
##              terminal nodes (tips)
##
########################################################
abstract type AbstractTip <: AbstractNode end

mutable struct ExtantTip <: AbstractTip
    index::Int64
    inbounds::Branch
    label::String
    sampling_probability::Float64
    is_fossil::Bool

    ExtantTip() = new()
end

mutable struct FossilTip <: AbstractTip
    index::Int64
    inbounds::Branch
    label::String

    FossilTip() = new()
end

mutable struct ExtinctionEvent <: AbstractTip
    index::Int64
    inbounds::Branch
    label::String

    ExtinctionEvent() = new()
end

########################################################
##
##              internal nodes
##
########################################################
abstract type InternalNode <: AbstractNode end
abstract type BranchingEvent <: InternalNode end

mutable struct Node <: BranchingEvent
    index::Int64
    inbounds::Branch
    children::Vector{Branch}
    Node() = new()
end

mutable struct Root <: BranchingEvent
    index::Int64
    children::Vector{Branch}
    Root() = new()
end

mutable struct SampledAncestor <: InternalNode
    index::Int64
    inbounds::Branch
    child::Branch

    SampledAncestor() = new()
end
