export writenewick
export newick

## write a newick file
@doc raw"""
writenewick(filename, data, rates)

writes a newick file with the rate values as comments

Example:
```julia
using Pesto

phy = readtree(Pesto.path("primates.tre"))
sampling_probability = 0.635
primates = SSEdata(phy, sampling_probability)

λ = [0.1, 0.2, 0.3, 0.4, 0.20]
μ = [0.05, 0.15, 0.05, 0.15, 0.25]
η = 1 / tree_length

model = BDSconstant(λ, μ, η)

rates = birth_death_shift(model, primates)

writenewick("/tmp/newick.tre", primates, rates)
```
"""
function writenewick(filename::String, data::SSEdata, rates::DataFrames.DataFrame)
    newick_string = newick(data, rates)

    open(filename, "w") do io
        write(io, newick_string)
        write(io, "\n")
    end
end

function node_data(rates::DataFrames.DataFrame)
    DataFrames.sort!(rates, :node)
    rates_plain = DataFrames.select(rates, DataFrames.Not(:node))
    
    res = []
    for row in eachrow(rates_plain)
        entries = String[]
        for (i, name) in enumerate(names(row))
            e = name * "=" * string(row[i])
            append!(entries, [e])
        end

        nd = String[]
        append!(nd, ["[&"])
        for (i, entry) in enumerate(entries)
            if i > 1
                append!(nd, [","])
            end
            append!(nd, [entry])
        end
        append!(nd, ["]"])

        append!(res, [*(nd...)])
    end
    return(res)
end

## create a newick string from the data object
## translated from R-package treeio: https://github.com/YuLab-SMU/treeio/blob/master/R/write-beast.R
function newick(data::SSEdata, rates::DataFrames.DataFrame)
    ancestors = make_ancestors(data)

    nd = node_data(rates)

    ind1 = [ancestors[i] for i in 1:ntip(data)]
    ind2 = [ancestors[i] for i in (ntip(data)+2):(nnode(data))]

    ind = vcat(ind1, 0, ind2)

    kids = make_descendants_nodes(data)

    root = getRoot(data.edges)
    desc = kids[root]
    s = String[]
    append!(s, ["("])
    n = ntip(data)

    for j in desc
        if j > n
            addinternal!(s, kids, nd, data, ind, j)
        else
            addterminal!(s, data, nd, ind[j])
        end

        if j != desc[length(desc)]
            append!(s, [","])
        end
    end
    append!(s, ["):0.0;"])
    
    newick = join(s)
    return(newick)
end

function addinternal!(s, kids, nd, data, ind, i)
    append!(s, ["("])

    desc = kids[i]

    for j in desc
        if j in data.edges[:,1]
            addinternal!(s, kids, nd, data, ind, j)
        else
            addterminal!(s, data, nd, ind[j])
        end

        if j != desc[length(desc)]
            append!(s, [","])
        end
    end

    append!(s, [")"])
    append!(s, [nd[i]])
    append!(s, [":"])
    append!(s, [string(data.branch_lengths[ind[i]])])
end

function addterminal!(s, data, nd, i)
    ii = data.edges[i,2]
    tl = data.tiplab[ii]

    append!(s, [tl])
    append!(s, [nd[ii]])
    append!(s, [":"])
    append!(s, [string(data.branch_lengths[i])])
end

function nnode(data)
    return (size(data.edges)[1]+1)
end

function ntip(data)
    return(length(data.tiplab))
end

function getRoot(edges)
    descendants = Set(edges[:,2])

    for node in edges[:,1]
        if node ∉ descendants
            return(node)
        end
    end
    throw("root not found")
end

