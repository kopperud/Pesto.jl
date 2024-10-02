## construct a pointer based tree from an APE style tree

function construct_tree(phy::phylo)
    
    tip_counter = Int64[0]

    root_index = length(phy.tip_label) + 1

    node_counter = Int64[root_index]

    
    root = Root()


    left_branch = construct_tree_po(phy, root, left_branch_index, tip_counter, node_counter) 

    
end

function construct_tree_po(
        phy,
        root::T,
        branch_index,
        tip_counter,
        node_counter
    ) where {T <: AbstractNode}


end
