export readtree, make_SSEdata, make_quantiles

function descendant_nodes(node, data)
    desc_edge_idxs = findall(data.edges[:,1] .== node)
    desc = data.edges[desc_edge_idxs,:]
    res = desc[:,2]
end

function parental_node(node, data)
    parental_edge_idx = findall(data.edges[:,2] .== node)
    parent_node = data.edges[parental_edge_idx,1][1]
    return(parent_node)
end

function make_quantiles(d, k)
    quantiles = zeros(k)
    for i in 1:k
        p = (i-0.5)/k
        quantiles[i] = Distributions.quantile(d, p)
    end
    return(quantiles)
end

function readtree(treefile)
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

    po <- postorder(phy)
    phy$po <- po
    """
    RCall.@rget phy
    return(phy)
end

function make_SSEdata(phy, datafile, ρ; include_traits = true)
   
    if include_traits
        df = CSV.File(datafile)
        trait_data = Dict(taxon => string(state) for (taxon, state) in zip(df[:Taxon], df[:state]))
    else
        trait_data = Dict(taxon => "?" for taxon in phy[:tip_label])
    end
    node_depth = phy[:node_depths]
    tiplab = phy[:tip_label]
    branching_times = phy[:branching_times]

    state_space = sort(unique(values(trait_data)))
    edges = convert.(Int64, phy[:edge])
    el = phy[:edge_length]
    po = phy[:po]

    data = SSEdata(state_space, trait_data, edges, tiplab, node_depth, ρ, el, branching_times, po)
    return(data)
end
