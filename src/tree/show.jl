function Base.Multimedia.display(root::Root)
    println("A phylogenetic root node (i.e. a tree) with $(length(root.children)) children.")
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
    println("$(length(branches)) branch structs.") 
end
