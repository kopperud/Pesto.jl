function find_one_tip(node::T) where {T <: InternalNode}
    find_one_tip(node.children[1].outbounds)
end

function find_one_tip(node::Tip)
    return(node)
end
