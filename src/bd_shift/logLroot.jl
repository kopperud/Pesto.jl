export logL_root

function number_of_states(model::TimevaryingModel)
    x = model.λ(0.0)
    n = length(x)
    return(n)
end

function number_of_states(model::ConstantModel)
    n = length(model.λ)
    return(n)
end

function get_speciation_rates(model::BDSconstant, t::Float64)
    return(model.λ)
end

function get_speciation_rates(model::BDStimevarying, t::Float64)
    return(model.λ(t)) 
end

function get_speciation_rates(model::FBDSconstant, t::Float64)
    return(model.λ)
end


function get_fossilization_rate(model::FBDSconstant, time::Float64)
    return(model.ψ)
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

    K = number_of_states(model)
    freqs = repeat([1.0 / K], K) ## prior on the root state

    # we condition the likelihood by
    #
    # * that there was a speciation event at the MRCA
    # * that the two lineages subtending from the MRCA 
    #        must have survived until the present
    λroot = get_speciation_rates(model, root_age)
    nonextinct = (1.0 .- E(root_age)).^2
    condition = λroot .* nonextinct

    D = D ./ condition
    prob = sum(freqs .* D)
    logL = log(prob) + sum(sf)
    return(logL)
end

function logL_root(model::Model, tree::Root)
    E = extinction_probability(model, tree)
    D, sf = postorder_async(model, tree, E)

    root_age = treeheight(tree)
    K = number_of_states(model)
    freqs = repeat([1.0 / K], K) ## prior on the root state

    # we condition the likelihood by
    #
    # * that there was a speciation event at the MRCA
    # * that the two lineages subtending from the MRCA 
    #        must have survived until the present
    λroot = get_speciation_rates(model, root_age)
    nonextinct = (1.0 .- E(root_age)).^2
    condition = λroot .* nonextinct

    D = D ./ condition
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
