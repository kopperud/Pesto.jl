
function contains_polytomies(phy::phylo)
    edge = phy[:edge]

    id = zeros(Int64, maximum(edge))
    for row in eachrow(edge)
        anc, dec = row
        id[anc] += 1
        id[dec] += 1
    end
    return(any(id .> 3))
end