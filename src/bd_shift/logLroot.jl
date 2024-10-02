export logL_root

function number_of_states(model::BDStimevarying)
    x = model.λ(0.0)
    n = length(x)
    return(n)
end

function number_of_states(model::BDSconstant)
    n = length(model.λ)
    return(n)
end

function get_speciation_rates(model::BDSconstant, t::Float64)
    return(model.λ)
end

function get_speciation_rates(model::BDStimevarying, t::Float64)
    return(model.λ(t)) 
end


function logL_root(model::Model, data::SSEdata; multithread = true)
    E = extinction_probability(model, data)

    if multithread
        D, sf = postorder_async(model, data, E)
    else
        D, sf = postorder_sync(model, data, E)
    end

    root_index = length(data.tiplab)+1
    root_age = data.node_depth[root_index]

    left_edge, right_edge = findall(data.edges[:,1] .== root_index)
#    λroot = get_speciation_rates(model, root_age)

    n = number_of_states(model)
    freqs = repeat([1.0 / n], n)

    # we divide by this to condition the probability density 
    # on that in order to have a tree in the first place, at 
    # least two lineages must have survived to the present.
    nonextinct = (1.0 .- E(root_age)).^2

    D = D ./ nonextinct
    prob = sum(freqs .* D)
    logL = log(prob) + sum(sf)
    return(logL)
end

function sselp(η, λ, μ, data)
    model = BDSconstant(λ, μ, η)

    logL_root(model, data)
end

function sselp_tv(η, λ, μ, data)
    model = BDStimevarying(λ, μ, t -> η)

    logL_root(model, data)
end
