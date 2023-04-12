## write a newick file
function write_tree(data, filename)
    newick_string = newick(data)

    touch(filename)
    io = open(filename, "w")

    write(io, newick_string)
    close(io)
end

function node_data(rates)
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
function newick(data, rates)
    ancestors = make_ancestors(data)

    nd = node_data(rates)

    ind1 = [ancestors[i] for i in 1:ntip(data)]
    ind2 = [ancestors[i] for i in (ntip(data)+2):(nnode(data))]

    ind = vcat(ind1, 0, ind2)

    kids = make_descendants_nodes(data)

    root = getRoot(data.edges)
    desc = kids[root]
    s = []
    append!(s, "(")
    n = ntip(data)

    for j in desc
        if j > n
            addinternal!(s, kids, nd, data, ind, j)
        else
            addterminal!(s, data, nd, ind[j])
        end

        if j != desc[length(desc)]
            append!(s, ",")
        end
    end
    append!(s, ");")
    
    newick = *(s...)
    return(newick)
end

function addinternal!(s, kids, nd, data, ind, i)
    append!(s, "(")

    desc = kids[i]

    for j in desc
        if j in data.edges[:,1]
            addinternal!(s, kids, nd, data, ind, j)
        else
            addterminal!(s, data, nd, ind[j])
        end

        if j != desc[length(desc)]
            append!(s, ",")
        end
    end

    append!(s, ")")
    append!(s, nd[ind[i]])
    append!(s, ":")
    append!(s, string(data.branch_lengths[ind[i]]))
end

function addterminal!(s, data, nd, i)
    ii = data.edges[i,2]

    tl = data.tiplab[ii]

    append!(s, tl)
    append!(s, ":")
    append!(s, string(data.branch_lengths[i]))
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

