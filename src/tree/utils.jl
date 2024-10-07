function find_one_tip(node::T) where {T <: InternalNode}
    find_one_tip(node.children[1].outbounds)
end

function find_one_tip(node::Tip)
    return(node)
end

function tip_labels(node::Root)
    labels = String[]
    tip_labels!(node, labels) 
    return(labels)
end

function tip_labels!(node::T, labels::Vector{String}) where {T <: InternalNode}
    for branch in node.children
        tip_labels!(branch.outbounds, labels)
    end
end

function tip_labels!(node::Tip, labels::Vector{String})
    push!(labels, node.label)
end

export get_branches

function get_branches(root::Root)
    branches = Branch[];

    get_branches!(root, branches)
    return(branches)
end

function get_branches!(node::T, branches) where {T <: InternalNode}
    for branch in node.children
        push!(branches, branch)
        get_branches!(branch.outbounds, branches)
    end
end

function get_branches!(node::Tip, branches)
    nothing
end

export number_of_branches

function number_of_branches(root::Root)
    n = [0]

    number_of_branches!(root, n)
    return(n[1])
end

function number_of_branches!(node::T, n::Vector{Int64}) where {T <: InternalNode}
    for branch in node.children
        n[1] += 1
        number_of_branches!(branch.outbounds, n)
    end
end

function number_of_branches!(tip::Tip, n::Vector{Int64})
    nothing
end

