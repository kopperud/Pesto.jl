## write tree to newick file

function writenewick(data)



end

function addinternal!(s, kids, i)
    append!(s, "(")

    desc = kids[i]
end

function addterminal!(s, i)

end

function nnode(data)
    return (size(data.edges)+1)
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


function writetree(data)
    s = String[]

    ancestors = data.edges[:,1]
    descendants = data.edges[:,2]

    #kids = zeros(Int64, size(data.edges))
    kids = make_descendants_nodes(data)

    n = length(data.tiplab)
    root = n+1

    desc_root = kids[root]

    append!(s, "(")
    for j in desc_root
        if j > n
            addinternal!(s, kids, i)
        else
            addterminal!(s, ind[j])
        end

        if j != desc[length(desc)]
            append!(s, ",")
        end
    end

    
end