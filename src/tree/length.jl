export treelength

function treelength(tree::Root)

    branch_lengths = Float64[]
    time = 0.0

    lengths!(tree, branch_lengths)

    tl = sum(branch_lengths)
    return(tl)
end

function lengths!(node::T, branch_lengths::Vector{Float64}) where {T <: InternalNode}
    for branch in node.children
        child_node = branch.outbounds
        push!(branch_lengths, branch.time)

        lengths!(child_node, branch_lengths)
    end
end

function lengths!(node::Tip, branch_lengths::Vector{Float64})
    nothing
end
