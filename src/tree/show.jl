function Base.Multimedia.display(tree::Root)
    println("A phylogenetic root node (i.e. a tree) with two branches, left and right.")
end

function Base.Multimedia.display(node::Node)
    println("An internal node, with three branches, inbounds, left and right.")
end

function Base.Multimedia.display(tip::Tip)
    println("A tip node with one incoming branch, representing \"$(tip.species_name)\".") 
end

function Base.Multimedia.display(branch::Branch)
    println("A branch with a parent node and a descendant node. The branch length is $(branch.times).")
end
