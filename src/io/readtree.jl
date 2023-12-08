export readtree

@doc raw"""
    readtree("/path/to/phylo.tre")

reads a NEXUS or Newick file using RCall and the R-package `ape`

Example:
```julia
using Pesto
phy = readtree(Pesto.path("primates.tre"))

display(phy)
```
"""
function readtree(treefile::String)
    s = readuntil(treefile, "(")
    isnexus = contains(s, "#NEXUS")

    RCall.@rput isnexus
    RCall.@rput treefile
    
    RCall.R"""
    library(ape)

    if(isnexus){
        phy <- read.nexus(treefile)
    }else{
        phy <- read.tree(treefile)
    }
    nde <- node.depth.edgelength(phy)
    node_depths <- max(nde) - nde
    phy$node_depths <- node_depths
    phy$branching_times <- branching.times(phy)
    phy$tip_label <- phy$tip.label

    po <- postorder(phy)
    phy$po <- po
    class(phy) <- "list"
    """
    RCall.@rget phy

    if !(phy[:branching_times] isa Vector)
        phy[:branching_times] = Float64[phy[:branching_times]]
    end

    r = phylo(
        phy[:edge],
        phy[:edge_length],
        phy[:Nnode],
        phy[:tip_label],
        phy[:node_depths],
        phy[:branching_times],
        phy[:po]
    )
    return(r)
end