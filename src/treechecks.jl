
function contains_polytomies(phy::phylo)
    edge = phy.edge

    id = zeros(Int64, maximum(edge))
    for row in eachrow(edge)
        anc, dec = row
        id[anc] += 1
        id[dec] += 1
    end
    return(any(id .> 3))
end

function is_ultrametric(phy::phylo)
    n_tips = length(phy.tip_label)

    @assert n_tips > 0 

    tip_depths = phy.node_depths[1:n_tips]

    ## because branch lengths are in floats
    ## the times can not be checked exactly 
    ## if they are equal to 0. Therefore some
    ## threshold is required
    threshold = 1e-4

    res = !any(abs.(tip_depths) .> threshold)
    return(res) 
end


