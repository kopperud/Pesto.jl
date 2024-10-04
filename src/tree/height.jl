export treeheight

function treeheight(tree::Root)

    node_depths = Float64[]
    time = 0.0

    depths!(node_depths, tree, time)

    height = maximum(node_depths)

    return(height)
end

function depths!(nd::Vector{Float64}, node::T, time::Float64) where {T <: InternalNode}
    for branch in node.children
        child_node = branch.outbounds      
        depths!(nd, child_node, time + branch.time)
    end
end

function depths!(nd::Vector{Float64}, node::Tip, time::Float64) 
    push!(nd, time)
end
