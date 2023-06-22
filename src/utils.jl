export readtree, make_SSEdata, make_quantiles, make_SSEdata2, allpairwise

function descendant_nodes(node, data)
    desc_edge_idxs = findall(data.edges[:,1] .== node)
    desc = data.edges[desc_edge_idxs,:]
    res = desc[:,2]
end

"""
    make_ancestors(data)

takes node indices as input, and returns edge indices
"""
function make_ancestors(data)
    ntip = length(data.tiplab)
    rootnode = ntip + 1
    maxnode = maximum(data.edges)

    ancestors = Dict(node => 0 for node in 1:maxnode if node != rootnode)

    for (i, row) in enumerate(eachrow(data.edges))
        anc, dec = row

        ancestors[dec] = i
    end
    return(ancestors)
end

"""
    make_descendants(data)

takes node indices as input, and returns edge indices
"""
function make_descendants(data)
    ntip = length(data.tiplab)
    rootnode = ntip + 1
    maxnode = maximum(data.edges)

    descendants = Dict(node => [] for node in rootnode:maxnode)

    for (i, row) in enumerate(eachrow(data.edges))
        anc, dec = row
        if anc > ntip
            append!(descendants[anc], i)
        end
    end
    return(descendants)
end

"""
    make_descendants(data)

takes node indices as input, and returns node indices
"""
function make_descendants_nodes(data)
    ntip = length(data.tiplab)
    rootnode = ntip + 1
    maxnode = maximum(data.edges)

    descendants = Dict(node => [] for node in rootnode:maxnode)

    for (i, row) in enumerate(eachrow(data.edges))
        anc, dec = row
        if anc > ntip
            append!(descendants[anc], dec)
        end
    end
    return(descendants)
end

function parental_node(node, data)
    parental_edge_idx = findall(data.edges[:,2] .== node)
    parent_node = data.edges[parental_edge_idx,1][1]
    return(parent_node)
end

function make_quantiles(d, k)
    quantiles = zeros(k)
    step = 0.5
    for i in 1:k
        p = (i-step)/k
        quantiles[i] = Distributions.quantile(d, p)
    end
    return(quantiles)
end

function make_quantiles2(d, k)
    ps = [(i-0.5)/k for i in 1:k]
    quantiles = Distributions.quantile.(d, ps)
    return(quantiles)
end

function make_quantiles3(d, k)
    ps = collect(range(0.0, 1.0; length = k+2))
    ps = ps[2:(length(ps)-1)]
    quantiles = Distributions.quantile.(d, ps)
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
    root_edge = phy[:root_edge]

    r = phylo(
        phy[:edge],
        phy[:edge_length],
        phy[:Nnode],
        phy[:tip_label],
        root_edge,
        phy[:node_depths],
        phy[:branching_times],
        phy[:po]
    )
    return(r)
end

function SSEdata(phy, ρ)
    make_SSEdata(phy, ρ)
end

function make_SSEdata(phy, datafile, ρ; include_traits = true)
    if contains_polytomies(phy)
        throw("Your tree is not a binary tree (it has hard polytomies). This program does not support trees with hard polytomies.")
    end
   
    if include_traits
        df = CSV.File(datafile)
        trait_data = Dict(taxon => string(state) for (taxon, state) in zip(df[:Taxon], df[:state]))
    else
        trait_data = Dict(taxon => "?" for taxon in phy.tip_label)
    end
    node_depth = phy.node_depths
    tiplab = phy.tip_label
    branching_times = phy.branching_times

    state_space = sort(unique(values(trait_data)))
    edges = convert.(Int64, phy.edge)
    el = phy.edge_length
    po = phy.po
    Nnode = phy.Nnode

    data = SSEdata(state_space, trait_data, edges, tiplab, node_depth, ρ, el, branching_times, po, Nnode)
    return(data)
end

@doc raw"""
    make_SSEdata(phy, ρ)

Example:
```julia
using Pesto
phy = readtree(Pesto.path("primates.tre")) 
ρ = 0.635  
data = make_SSEdata(phy, ρ) 
```
"""
function make_SSEdata(phy, ρ)
    trait_data = Dict(taxon => "?" for taxon in phy.tip_label)

    node_depth = phy.node_depths
    tiplab = phy.tip_label
    branching_times = phy.branching_times

    state_space = NaN
    edges = convert.(Int64, phy.edge)
    el = phy.edge_length
    po = phy.po
    if any(el .< 0)
        throw(error("Tree includes negative branch lengths."))
    end
    Nnode = phy.Nnode

    data = SSEdata(state_space, trait_data, edges, tiplab, node_depth, ρ, el, branching_times, po, Nnode)
    return(data)
end

function partition_postorder_indices(data)
    ancestor_node = Dict(val => key for (key, val) in eachrow(data.edges))

    d = Dict(node => 0 for node in data.edges[:,2])
    parents = collect(1:length(data.tiplab))
    res = [parents]

    root_node = length(data.tiplab)+1
    for i in 1:maximum(data.edges)
        for node in parents
            parent = ancestor_node[node]
            if parent != root_node
                d[parent] += 1
            end 
        end

        # find the twos
        parents = Int64[]
        for (key, val) in d
            if val == 2
                append!(parents, key)
                delete!(d, key)
            end
        end

        if isempty(parents)
            break
        end
        append!(res, [parents])
    end

    return(res)
end


@doc raw"""
    allpairwise(λ, μ)

Example:
```julia
using Pesto

lambda = [0.2, 0.3]
mu = [0.05, 0.10, 0.15, 0.20] 

λ, μ = allpairwise(lambda, mu)
```
"""
function allpairwise(xs, ys)
    ny = length(xs)
    nx = length(ys)

    k = ny * nx

    λ = zeros(eltype(xs), k)
    μ = zeros(eltype(ys), k)
    
    for (i, (x, y)) in enumerate(Iterators.product(xs, ys))
        λ[i] = x
        μ[i] = y
    end

    return(λ, μ)
end

#function (from = 1, to = 1e+05, length.out = 6) 
#    {
#        exp(seq(log(from), log(to), length.out = length.out))
#    }

function lrange(from, to; length = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end

