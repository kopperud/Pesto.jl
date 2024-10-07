function Base.Multimedia.display(root::Root)
    tl = tip_labels(root)
    too_long = length(tl) > 10

    if too_long 
        tl = tl[1:10]
        #tl_short = tl_short[1:50]
    end

    tl_short = join(tl, ", ")

    if too_long
        tl_short = tl_short * "..."
    end

    height = round(treeheight(root); digits = 2)

    println("A phylogenetic root node (i.e. a tree) with $(length(root.children)) children, height $height and tip labels $tl_short")
end

function Base.Multimedia.display(node::Node)
    println("An internal node, with an incoming branch, and $(length(node.children)) children.")
end

function Base.Multimedia.display(tip::Tip)
    println("A tip node with one incoming branch, representing \"$(tip.label)\", with sampling probability $(tip.sampling_probability).") 
end

function Base.Multimedia.display(branch::Branch)
    println("A branch with a parent node and a descendant node. The branch length is $(branch.time).")
end

function Base.Multimedia.display(branches::Vector{Branch})
    println("A vector of $(length(branches)) branches.") 
end
