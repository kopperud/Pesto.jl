export logL_root


#=
function number_of_states(model::TimevaryingModel)
    x = model.λ(0.0)
    n = length(x)
    return(n)
end
=#

function number_of_states(model::MultiStateModel)
    n = length(model.λ)
    return(n)
end

function number_of_states(model::BhtvDhtvModel)
    n = length(model.λ(0.0))
    return(n)
end

function number_of_states(model::UniStateModel)
    return(1)
end

#=
function get_speciation_rates(model::BDSconstantQ, t::Float64)
    return(model.λ)
end
=#

#=
function get_speciation_rates(model::BDStimevarying, t::Float64)
    return(model.λ(t)) 
end
=#


function logL_root(model::Model, data::SSEdata; condition = [:survival, :mrca])
    u, sf = postorder_async(model, data)

    E = u[:,1]
    D = u[:,2]

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
    if :survival in condition
        nonextinct = (1.0 .- E).^2
        D = D ./ nonextinct
    end

    if :mrca in condition
        λroot = get_speciation_rates(model, root_age)
        D = D ./ λroot
    end

    prob = sum(freqs .* D)
    logL = log(prob) + sum(sf)
    return(logL)
end

function logL_root(model::Model, tree::Root; condition = [:survival, :mrca])
    u, sf = postorder_async(model, tree)

    E, D = eachcol(u)

    root_age = treeheight(tree)
    K = number_of_states(model)

    freqs = repeat([1.0 / K], K) ## prior on the root state

    # we condition the likelihood by
    #
    # * that there was a speciation event at the MRCA
    # * that the two lineages subtending from the MRCA 
    #        must have survived until the present
    if :survival in condition
        nonextinct = (1.0 .- E).^2
        D = D ./ nonextinct
    end

    if :mrca in condition
        λroot = get_speciation_rates(model, root_age)
        D = D ./ λroot
    end

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
