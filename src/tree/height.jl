export treeheight
export node_depths

function treeheight(tree::Root)

    nd = node_depths(tree)
    height = maximum(nd)

    return(height)
end

function node_depths(tree::Root)
    nd = Float64[]

    time = 0.0
    depths!(nd, tree, time)

    return(nd)
end

function depths!(
        nd::Vector{Float64}, 
        node::T, 
        time::Float64
    ) where {T <: BranchingEvent}

    push!(nd, time)

    for branch in node.children
        child_node = branch.outbounds      
        depths!(nd, child_node, time + branch.time)
    end
end

function depths!(nd::Vector{Float64}, node::T, time::Float64) where {T <: AbstractTip}
    push!(nd, time)
end

###
function branch_times(tree::Root)
    t_old = treeheight(tree)

    n_branches = number_of_branches(tree)

    ## t_old and t_young for all branches
    times = zeros(Float64, n_branches, 2)

    branch_times!(tree, times, t_old)

    return(times)
end

function branch_times!(node::T, times::Matrix{Float64}, t_old::Float64) where {T <: BranchingEvent}
    for branch in node.children 
        t_young = t_old - branch.time
        times[branch.index,1] = t_old
        times[branch.index,2] = t_young

        branch_times!(branch.outbounds, times, t_young)
    end
end

function branch_times!(tip::T, times::Matrix{Float64}, t_old::Float64) where {T <: AbstractTip}
    nothing
end






