function find_one_extant_tip(root::Root)
    tips = ExtantTip[]

    find_one_extant_tip!(root, tips)

    return(tips[1])
end

function find_one_extant_tip!(
        node::T, 
        tips::Vector{ExtantTip}
    ) where {T <: BranchingEvent}

    for branch in node.children
        if length(tips) == 0
            find_one_extant_tip!(branch.outbounds, tips)
        end
    end
end

function find_one_extant_tip!(node::ExtantTip, tips::Vector{ExtantTip})
    push!(tips, node)
end

function find_one_extant_tip!(node::FossilTip, tips::Vector{ExtantTip})
    nothing
end


function tip_labels(node::Root)
    labels = String[]
    tip_labels!(node, labels) 
    return(labels)
end

function tip_labels!(node::T, labels::Vector{String}) where {T <: BranchingEvent}
    for branch in node.children
        tip_labels!(branch.outbounds, labels)
    end
end

function tip_labels!(node::T, labels::Vector{String}) where {T <: AbstractTip}
    push!(labels, node.label)
end

export get_branches

function get_branches(root::Root)
    branches = Branch[];

    get_branches!(root, branches)
    return(branches)
end

function get_branches!(node::T, branches) where {T <: BranchingEvent}
    for branch in node.children
        push!(branches, branch)
        get_branches!(branch.outbounds, branches)
    end
end

function get_branches!(node::T, branches) where {T <: AbstractTip}
    nothing
end

export number_of_branches

function number_of_branches(root::Root)
    n = [0]

    number_of_branches!(root, n)
    return(n[1])
end

function number_of_branches!(node::T, n::Vector{Int64}) where {T <: BranchingEvent}
    for branch in node.children
        n[1] += 1
        number_of_branches!(branch.outbounds, n)
    end
end

function number_of_branches!(tip::T, n::Vector{Int64}) where {T <: AbstractTip}
    nothing
end

