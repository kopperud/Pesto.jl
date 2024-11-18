## construct a pointer based tree from an APE style tree

export construct_tree

function construct_tree(phy::phylo, sampling_probability::Float64)
    root_index = length(phy.tip_label) + 1

    root = Root()
    root.children = Vector{Branch}()
    root.index = root_index

    descendant_edges = findall(phy.edge[:,1] .== root_index)

    @assert length(descendant_edges) == 2

    left_branch_index, right_branch_index = descendant_edges
    left_node_index  = phy.edge[descendant_edges[1],2]
    right_node_index = phy.edge[descendant_edges[2],2]

    if left_node_index > root_index
        left_branch = construct_tree_po(phy, root, left_branch_index, sampling_probability)
    else
        left_branch = construct_tree_tip(phy, root, left_branch_index, sampling_probability)
    end

    if right_node_index > root_index
        right_branch = construct_tree_po(phy, root, right_branch_index, sampling_probability)
    else
        right_branch = construct_tree_tip(phy, root, right_branch_index, sampling_probability)
    end
    
    return(root)
end

function construct_tree_po(
        phy::phylo,
        parent_node::T,
        branch_index::Int64,
        sampling_probability::Float64,
    ) where {T <: AbstractNode}


    branch = Branch()
    branch.time = phy.edge_length[branch_index]
    branch.inbounds = parent_node
    branch.index = branch_index 
    push!(parent_node.children, branch)

    node = Node()
    node.inbounds = branch
    node.children = Vector{Branch}()
    branch.outbounds = node

    node_index = phy.edge[branch_index,2]
    node.index = node_index


    root_index = length(phy.tip_label) + 1
    descendant_edges = findall(phy.edge[:,1] .== node_index)
    @assert length(descendant_edges) == 2

    left_branch_index, right_branch_index = descendant_edges
    left_node_index  = phy.edge[descendant_edges[1],2]
    right_node_index = phy.edge[descendant_edges[2],2]

    if left_node_index > root_index
        left_branch = construct_tree_po(phy, node, left_branch_index, sampling_probability)
    else
        left_branch = construct_tree_tip(phy, node, left_branch_index, sampling_probability)
    end


    if right_node_index > root_index
        construct_tree_po(phy, node, right_branch_index, sampling_probability)
    else
        construct_tree_tip(phy, node, right_branch_index, sampling_probability)
    end

    nothing
end

function construct_tree_tip(phy, parent_node, branch_index, sampling_probability)
    branch = Branch()
    branch.time = phy.edge_length[branch_index]
    branch.inbounds = parent_node
    branch.index = branch_index
    push!(parent_node.children, branch)

    tip = ExtantTip()
    tip_index = phy.edge[branch_index,2]
    tip.index = tip_index
    tip.label = phy.tip_label[tip_index]
    tip.inbounds = branch
    tip.sampling_probability = sampling_probability
    branch.outbounds = tip

    nothing
end


function assign_fossils!(
        root::Root,
        threshold::Float64,
    )
    time = treeheight(root)

    assign_fossils!(root, threshold, time)
end

function assign_fossils!(
        node::T,
        threshold::Float64,
        time::Float64,
    ) where {T <: BranchingEvent}

    for branch in node.children  
        child_node = branch.outbounds
        assign_fossils!(child_node, threshold, time - branch.time)
    end
    println()
end

function assign_fossils!(
        tip::ExtantTip,
        threshold::Float64,
        time::Float64,
    )
    if !(time < threshold)
        
        ft = FossilTip()
        ft.index = tip.index
        ft.inbounds = tip.inbounds 
        ft.label = tip.label
        parent_branch = tip.inbounds
        parent_branch.outbounds = ft

        print(".")
    end
end


