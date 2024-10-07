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
