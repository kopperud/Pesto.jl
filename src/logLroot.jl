export logL_root

function logL_root(model::SSEconstant, data::SSEdata)
    E = extinction_probability(model, data)
    D_ends, sf = postorder_nosave(model, data, E)
    root_index = length(data.tiplab)+1
    root_age = data.node_depth[root_index]

    left_edge, right_edge = findall(data.edges[:,1] .== root_index)
    D_left = D_ends[left_edge,:]
    D_right = D_ends[right_edge,:]
    D = D_left .* D_right .* model.λ

    #    freqs = [0.5, 0.5]
    n = length(model.λ)
    freqs = repeat([1.0 / n], n)

    # we divide by this to condition the probability density 
    # on that in order to have a tree in the first place, at 
    # least two lineages must have survived to the present.
    nonextinct = (1.0 .- E(root_age)).^2
    #D = arr[end,:,2]
    # Why do we divide by (1-E)^2 * λ, and not just (1-E)^2 ?
#    D = D ./ (nonextinct .* [λ(root_age) for λ in model.λ])
    D = D ./ (nonextinct .* model.λ)
    prob = sum(freqs .* D)
    logL = log(prob) + sum(sf)
end
