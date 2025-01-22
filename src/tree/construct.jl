## construct a pointer based tree from an APE style tree

export construct_tree

function construct_tree(phy::phylo, sampling_probability::Float64; fossil_threshold = 0.5)
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

    ## assign fossils
    assign_fossil_tips!(root, fossil_threshold)
    ## assign fossils
    assign_sampled_ancestors!(root)
    
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


function assign_fossil_tips!(
        root::Root,
        threshold::Float64,
    )
    time = treeheight(root)

    n = assign_fossil_tips!(root, threshold, time)
    if n > 0
        println("Assigned $n tips as fossil samples.")
    end
end

function assign_fossil_tips!(
        node::T,
        threshold::Float64,
        time::Float64,
    ) where {T <: BranchingEvent}

    n = 0
    for branch in node.children  
        child_node = branch.outbounds
        n += assign_fossil_tips!(child_node, threshold, time - branch.time)
    end

    return(n)
end

function assign_fossil_tips!(
        tip::ExtantTip,
        threshold::Float64,
        time::Float64,
    )
    n = 0

    if !(time < threshold)
        ft = FossilTip()
        ft.index = tip.index
        ft.inbounds = tip.inbounds 
        ft.label = tip.label
        parent_branch = tip.inbounds
        parent_branch.outbounds = ft
        n += 1
    end
    return(n)
end

function assign_fossil_tips!(
        tip::FossilTip,
        threshold::Float64,
        time::Float64,
    )
    return(0)
end

########################
##
##   turn 0-edge branches in to fossil sampled ancestors
##
#########################
function get_sister_branch(branch::Branch)
    parent_node = branch.inbounds

    branches = parent_node.children

    is_other = [branch != b for b in branches]
    sister_branch = branches[findfirst(is_other)]

    return(sister_branch)
end

function assign_sampled_ancestors!(root::Root)
    n = 0
    for branch in root.children
        n = assign_sampled_ancestors!(branch)
    end

    if n > 0
        println("Assigned $n internal nodes as being sampled ancestors.")
    end
end


function assign_sampled_ancestors!(branch::Branch)
    threshold = 1e-8

    is_zero = abs(branch.time) < threshold

    n = 0
    if is_zero
        sister_branch = get_sister_branch(branch)

        parent_node = branch.inbounds
        if parent_node isa Root
            error("sampled ancestor descending directly from root node is not supported.")
        else
            parent_branch = parent_node.inbounds

            sampled_ancestor = SampledAncestor()
            sampled_ancestor.index = 9999
            sampled_ancestor.inbounds = parent_branch
            parent_branch.outbounds = sampled_ancestor

            sampled_ancestor.child = sister_branch
            sister_branch.inbounds = sampled_ancestor
            n += 1
        end
    end
    n += assign_sampled_ancestors!(branch.outbounds)
    return(n)
end

function assign_sampled_ancestors!(sampled_ancestor::SampledAncestor)
    n = assign_sampled_ancestors!(sampled_ancesor.child)
    return(n)
end

function assign_sampled_ancestors!(node::BranchingEvent)
    n = 0
    for branch in node.children
        n = assign_sampled_ancestors!(branch)
    end
    return(n)
end

function assign_sampled_ancestors!(tip::ExtantTip)
    n = 0
    return(n)
end
function assign_sampled_ancestors!(tip::FossilTip)
    n = 0
    return(n)
end



#=
function assign_sampled_ancestors!(tree::Root)
    threshold = 1e-8

    n = 0
    branches = tree.children
    @assert length(branches) == 2

    is_zero = [abs(branch.time) < threshold for branch in branches]
    if sum(is_zero) > 1
        error("more than two branches descending from the same node were zero-length. this does not make sense.")
    elseif sum(is_zero) == 1
        which_branch = findfirst(is_zero)
        sister_branch = findfirst(.!is_zero)
        child_node = branches[which_branch].outbounds

        if child_node isa ExtantTip
            sampled_ancestor = SampledAncestor() 
            sampled_ancestor.index = 9999
            sampled_ancestor.inbounds = sister_branch
            sister_branch.outbounds = sampled_ancestor
            #sampled_ancestor.outbounds = 
        else
            error("the child node of the zero-length branch is not a fossil tip. this does not make sense.")
        end
    end


    #if any([abs(branch.time) < threshold for branch in branches])
    for branch in tree.children
        n = assign_sampled_ancestors!(child_node)
    end

    if n > 0
        println("Assigned $n internal nodes as being sampled ancestors.")
    end
end
=#

