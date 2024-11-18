export make_SSEdata, make_quantiles, SSEdata, allpairwise, lrange, make_descendants, make_ancestors, make_descendants_nodes

function descendant_nodes(node, data)
    desc_edge_idxs = findall(data.edges[:,1] .== node)
    desc = data.edges[desc_edge_idxs,:]
    res = desc[:,2]
end

@doc raw"""
    make_ancestors(data)

takes node indices as input, and returns edge indices

Example:
```julia
using Pesto
phy = readtree(Pesto.path("primates.tre"))
sampling_probability = 0.635
data = make_SSEdata(phy, sampling_probability)

make_ancestors(data)
```
with result
```julia
Dict{Int64, Int64} with 464 entries:
  56  => 115
  35  => 71
  425 => 379
  ⋮   => ⋮
```
"""
function make_ancestors(data::SSEdata)
    ntip = length(data.tiplab)
    rootnode = ntip + 1
#    maxnode = maximum(data.edges)
    maxnode = data.Nnode + length(data.tiplab)    

    ancestors = Dict{Int64, Int64}()    
    #ancestors = zeros(Int64, maxnode)

    #for (i, row) in enumerate(eachrow(data.edges))
    for i in 1:size(data.edges)[1]
        dec = data.edges[i,2]

        ancestors[dec] = i
    end
    return(ancestors)
end

@doc raw"""
    make_descendants(data)

takes node indices as input, and returns edge indices
Example:
```julia
using Pesto
phy = readtree(Pesto.path("primates.tre"))
sampling_probability = 0.635
data = make_SSEdata(phy, sampling_probability)

make_descendants(data)
```
with result
```julia
Dict{Int64, Vector{Any}} with 232 entries:
  402 => [330, 331]
  413 => [357, 360]
  425 => [380, 381]
  ⋮   => ⋮
```
"""
function make_descendants(data::SSEdata)
    ntip = length(data.tiplab)
    rootnode = ntip + 1
    maxnode = maximum(data.edges)

    descendants = Dict(node => Int64[] for node in rootnode:maxnode)

    for (i, row) in enumerate(eachrow(data.edges))
        anc, dec = row
        if anc > ntip
            append!(descendants[anc], i)
        end
    end
    return(descendants)
end

function make_descendants(edges::Array{Int64,2})
    ntip = size(edges)[1]÷2 + 1
    rootnode = ntip + 1
    maxnode = maximum(edges)

    descendants = Dict(node => [] for node in rootnode:maxnode)

    for (i, row) in enumerate(eachrow(edges))
        anc, dec = row
        if anc > ntip
            append!(descendants[anc], i)
        end
    end
    return(descendants)
end

@doc raw"""
    make_descendants_nodes(data)

takes node indices as input, and returns node indices
Example:
```julia
using Pesto
phy = readtree(Pesto.path("primates.tre"))
sampling_probability = 0.635
data = make_SSEdata(phy, sampling_probability)

make_descendants_nodes(data)
```
"""
function make_descendants_nodes(data::SSEdata)
    ntip = length(data.tiplab)
    rootnode = ntip + 1
    maxnode = maximum(data.edges)

    #descendants = Dict(node => Int64[] for node in rootnode:maxnode)
    descendants = Dict{Int64, Vector{Int64}}()
    for node in rootnode:maxnode
        descendants[node] = Int64[]
    end


    for row in eachrow(data.edges)
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

@doc raw"""
    SSEdata(phy, sampling_probability)

Example:
```julia
using Pesto
phy = readtree(Pesto.path("primates.tre")) 
sampling_probability = 0.635  
data = SSEdata(phy, sampling_probability) 
```
"""
function SSEdata(phy::phylo, sampling_probability::Float64)
    make_SSEdata(phy, "", sampling_probability; include_traits = false)
end

function make_SSEdata(phy::phylo, datafile::String, sampling_probability::Float64; include_traits = true)
    if contains_polytomies(phy)
        #throw("Your tree is not a binary tree (it has hard polytomies). This program does not support trees with hard polytomies.")
        throw(PolytomyError())
    end

    if !is_ultrametric(phy)
        @warn "Your tree appears to not be ultrametric. Check if this is a rounding error issue, or if it really is not ultrametric. Pesto only works for ultrametric trees."
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

    if any(el .< 0)
        throw(NegativeBranchError())
    end

    if any(el .< 1e-6)
        @warn "Your tree has some branch lengths that are very short (<0.000001). Check if this is really is a good representation of the divergence times."
    end

    data = SSEdata(state_space, trait_data, edges, tiplab, node_depth, sampling_probability, el, branching_times, po, Nnode)
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
with result
```julia
([0.2, 0.3, 0.2, 0.3, 0.2, 0.3, 0.2, 0.3], [0.05, 0.05, 0.1, 0.1, 0.15, 0.15, 0.2, 0.2])
```
"""
function allpairwise(xs, ys)
    ny = length(xs)
    nx = length(ys)

    k = ny * nx

    λ = zeros(Base.eltype(xs), k)
    μ = zeros(Base.eltype(ys), k)
    
    for (i, (x, y)) in enumerate(Iterators.product(xs, ys))
        λ[i] = x
        μ[i] = y
    end

    return(λ, μ)
end

export alltriples

function alltriples(x1::Vector{T}, x2::Vector{T}, x3::Vector{T}) where {T <: Real}
    n1 = length(x1)
    n2 = length(x2)
    n3 = length(x3)

    k = n1 * n2 * n3 

    λ = zeros(T, k)
    μ = zeros(T, k)
    ψ = zeros(T, k)

    for (i, (x, y, z)) in enumerate(Iterators.product(x1, x2, x3))
        λ[i] = x
        μ[i] = y
        ψ[i] = z
    end

    return(λ, μ, ψ)
end

#function Qmatrix(λ::Vector{T}, μ::Vector{T}, ψ::Vector{T}, α, β, γ) where {T <: Real}
function Qmatrix(r::Vector{T}, ψ::Vector{T}, α::T, β::T) where {T <: Real}
    ## backtransform
    #λ = r ./ (1.0 - ϵ)
    #μ = λ .- r
    
    k = reduce(*, map(length, (r, ψ)))

    nr = length(r)
    nψ = length(ψ)

    small_α = α / (nr - 1)
    small_β = β / (nψ - 1)

    ## empty flat vectors
    A = zeros(T, k)
    B = zeros(T, k)
    #C = zeros(T, k)

    ## compute all pairs, populate in vectors
    for (i, (a, b)) in enumerate(Iterators.product(r, ψ))
        A[i] = a
        B[i] = b
        #C[i] = c
    end
   
    ## Q matrix
    Q = zeros(T, k, k)
    Qα = zeros(Int64, k, k)
    Qβ = zeros(Int64, k, k)
    #Qγ = zeros(Int64, k, k)
    
    for i in 1:k
        for j in 1:k
            
            ## only one of r or ψ is permitted to change
            ## otherwise is it not assigned (default to zero)
            n_changes = (A[i] != A[j]) + (B[i] != B[j])# + (C[i] != C[j])

            if n_changes == 1
                if A[i] != A[j]
                    Q[i,j] = small_α
                    Qα[i,j] = 1
                else
                #elseif B[i] != B[j]
                    Q[i,j] = small_β
                    Qβ[i,j] = 1
                end
                #else
                #    Q[i,j] = small_γ
                #    Qγ[i,j] = 1
                #end
            end
        end
    end

    for i in 1:k
        #Q[i,i] = -sum(Q[:,i])
        Q[i,i] = -(α+β)
    end

    Qs = SparseArrays.sparse(Q)
    Qαs = SparseArrays.sparse(Qα)
    Qβs = SparseArrays.sparse(Qβ)
    #Qγs = SparseArrays.sparse(Qγ)

    return(Qs, Qαs, Qβs) #, Qγs)
end

@doc raw"""
    lrange(from, to, length)

Similar to `range`, but with proportional spacing.

Example:
```julia
using Pesto

lrange(0.001, 100.0, 6)
```
with result
```julia
6-element Vector{Float64}:
   0.0010000000000000002
   0.010000000000000004
   0.10000000000000002
   1.0000000000000004
  10.000000000000002
 100.00000000000004
```
"""
function lrange(from::Float64, to::Float64, length::Int64 = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end

function getpar(x::Float64)
    return(x)
end

function getpar(x::ForwardDiff.Dual)
    return(getpar(x.value)) ## recursive incase of higher-order derivatives
end



notneg(u,p,t) = any(x->x<0,u)

function eltype(model::BDSconstant)
    return(typeof(model.η))
end
function eltype(model::BDSconstantQ)
    return(typeof(model.Q[1,1]))
end
function eltype(model::BDStimevarying)
    return(typeof(model.η(0.0)))
end
function eltype(model::FBDSconstant)
    return(typeof(model.α))
end

