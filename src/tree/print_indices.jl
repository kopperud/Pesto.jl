function print_indices(node::T) where {T <: AbstractNode}
    #println("this is a node with index $(node.index)")

    if !(node isa Tip)
        for branch in node.children
            print_indices(branch)
        end
    end
end

function print_indices(branch::Branch)
    println("this is a branch with index $(branch.index)")

    print_indices(branch.outbounds)
end


