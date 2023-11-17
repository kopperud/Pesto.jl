export logL_root

function number_of_states(model::SSEtimevarying)
    x = model.λ(0.0)
    n = length(x)
    return(n)
end

function number_of_states(model::SSEconstant)
    n = length(model.λ)
    return(n)
end

function get_speciation_rates(model::SSEconstant, t::Float64)
    return(model.λ) 
end

function get_speciation_rates(model::SSEtimevarying, t::Float64)
    return(model.λ(t)) 
end

function logL_root(model::SSE, data::SSEdata)
    E = extinction_probability(model, data)
    D_ends, sf = postorder_nosave(model, data, E)
    root_index = length(data.tiplab)+1
    root_age = data.node_depth[root_index]

    left_edge, right_edge = findall(data.edges[:,1] .== root_index)
    D_left = D_ends[left_edge,:]
    D_right = D_ends[right_edge,:]
#    λroot = get_speciation_rates(model, root_age)
    D = D_left .* D_right 

    n = number_of_states(model)
    freqs = repeat([1.0 / n], n)

    # we divide by this to condition the probability density 
    # on that in order to have a tree in the first place, at 
    # least two lineages must have survived to the present.
    nonextinct = (1.0 .- E(root_age)).^2

    # Why do we divide by (1-E)^2 * λ, and not just (1-E)^2 ?
    D = D ./ nonextinct
    prob = sum(freqs .* D)
    logL = log(prob) + sum(sf)
end

function sselp(η, λ, μ, data)
    model = SSEconstant(λ, μ, η)

    logL_root(model, data)
end

function sselp_tv(η, λ, μ, data)
    model = SSEtimevarying(λ, μ, t -> η)

    logL_root(model, data)
end
