export tree_rates
export ancestral_state_probabilities

#https://en.wikipedia.org/wiki/Gaussian_quadrature
function quadrature(f, t0, t1, x, w)
    ## change of interval from t0 => t1 to -1 => 1
    g(x) = ((t1 - t0) / 2.0) * f(((t1 - t0)/2.0)*x + (t1 + t0)/2.0)
    
    I = LinearAlgebra.dot(w, g.(x))
    return(I)
end

function meanbranch(f, t0, t1, x, w)
    # I is the integral
    I = quadrature(f, t0, t1, x, w)

    # divide by time interval, since we want the average rate
    res = I / (t1 - t0)
    return(res)
end

function tree_rates(data, model; n = 10)
    Ds, Fs = backwards_forwards_pass(model, data);
    Ss = ancestral_state_probabilities(data, Ds, Fs);

    tree_rates(data, model, Fs, Ss; n = n)
end

function tree_rates(data::SSEdata, model::T, Fs, Ss; n = 10) where {T <: MultiStateModel}
    #rates = zeros(size(data.edges)[1], 8)
    x, w = FastGaussQuadrature.gausslegendre(n)

    n_branches = number_of_branches(data) 
    keys = [
         "mean_lambda", "mean_mu", "mean_netdiv", "mean_relext",
         "delta_lambda", "delta_mu", "delta_netdiv", "delta_relext", 
        ]
    if hasproperty(model, :ψ)
        push!(keys, "mean_psi")
        push!(keys, "delta_psi")
    end

    res = Dict{String,Vector{Float64}}()
    for key in keys
        res[key] = zeros(n_branches)
    end
    
    Threads.@threads for i = 1:n_branches
        t0, t1 = extrema(Fs[i].t)
        ## t0 is youngest, t1 is oldest
        
        ## posterior mean rate, and
        ## difference from oldest to youngest point on branch
        
        ## speciation rate
        res["mean_lambda"][i] = meanbranch(t -> LinearAlgebra.dot(get_speciation_rates(model, t), Ss[i](t)), t0, t1, x, w)
        lambda_young = LinearAlgebra.dot(get_speciation_rates(model, t0), Ss[i](t0))
        lambda_old   = LinearAlgebra.dot(get_speciation_rates(model, t1), Ss[i](t1))
        res["delta_lambda"][i] = lambda_young - lambda_old
        ## extinction rate
        res["mean_mu"][i] = meanbranch(t -> LinearAlgebra.dot(get_extinction_rates(model, t), Ss[i](t)), t0, t1, x, w)
        mu_young = LinearAlgebra.dot(get_extinction_rates(model, t0), Ss[i](t0))
        mu_old   = LinearAlgebra.dot(get_extinction_rates(model, t1), Ss[i](t1))
        res["delta_mu"][i] = mu_young - mu_old

        ## net-diversification rate
        res["mean_netdiv"][i] = meanbranch(t -> LinearAlgebra.dot(get_speciation_rates(model, t) .- get_extinction_rates(model, t), Ss[i](t)), t0, t1, x, w)
        netdiv_young = LinearAlgebra.dot(get_speciation_rates(model, t0) .- get_extinction_rates(model, t0), Ss[i](t0))
        netdiv_old   = LinearAlgebra.dot(get_speciation_rates(model, t1) .- get_extinction_rates(model, t1), Ss[i](t1))
        res["delta_netdiv"][i] = netdiv_young - netdiv_old 

        ## relative extinction rate (μ/λ)
        res["mean_relext"][i] = meanbranch(t -> LinearAlgebra.dot(get_extinction_rates(model, t) ./ get_speciation_rates(model, t), Ss[i](t)), t0, t1, x, w)
        relext_young = LinearAlgebra.dot(get_extinction_rates(model, t0) ./ get_speciation_rates(model, t0), Ss[i](t0)) 
        relext_old   = LinearAlgebra.dot(get_extinction_rates(model, t1) ./ get_speciation_rates(model, t1), Ss[i](t1)) 
        res["delta_relext"][i] = relext_young - relext_old

        ## only if the model is an FBD model
        if hasproperty(model, :ψ)
            res["mean_psi"][i] = meanbranch(t -> LinearAlgebra.dot(model.ψ, Ss[i](t)), t0, t1, x, w)
            res["delta_psi"][i] = LinearAlgebra.dot(model.ψ, Ss[i](t0)) - LinearAlgebra.dot(model.ψ, Ss[i](t1))
        end
    end

    #node = data.edges[:,2]
    edge = 1:n_branches

    df = DataFrames.DataFrame(res)
    df[!, "node"] = get_node_indices(data)
    df[!, "edge"] = edge

    tiplabs = tip_labels(data)
    root_index = length(tiplabs)+1

    root_row = Dict(key => NaN for key in keys)
    root_row["node"] = root_index
    root_row["edge"] = 0

    df = push!(df, root_row)
    #push!(df, vcat(repeat([NaN], length(names)), root_index, 0))
    return(df)
end

function tree_rates(tree::Root, model::T, Fs, Ss; n = 10) where {T <: MultiStateModel}
    #rates = zeros(size(data.edges)[1], 8)
    x, w = FastGaussQuadrature.gausslegendre(n)

    n_branches = number_of_branches(tree)
    keys = [
         "mean_lambda", "mean_mu", "mean_netdiv", "mean_relext",
         "delta_lambda", "delta_mu", "delta_netdiv", "delta_relext", 
        ]
    if hasproperty(model, :ψ)
        push!(keys, "mean_psi")
        push!(keys, "delta_psi")
    end

    branches = get_branches(tree)
    maxindex = maximum(branch.index for branch in branches)

    res = Dict{String,Vector{Float64}}()
    for key in keys
        res[key] = zeros(maxindex)
    end
    
    #Threads.@threads for i = 1:n_branches

    #for (i, _) in branches
    for branch in branches
        i = branch.index

        t0, t1 = extrema(Fs[i].t)
        ## t0 is youngest, t1 is oldest
        
        ## posterior mean rate, and
        ## difference from oldest to youngest point on branch
        
        ## speciation rate
        res["mean_lambda"][i] = meanbranch(t -> LinearAlgebra.dot(model.λ, Ss[i](t)), t0, t1, x, w)
        res["delta_lambda"][i] = LinearAlgebra.dot(model.λ, Ss[i](t0)) - LinearAlgebra.dot(model.λ, Ss[i](t1))

        ## extinction rate
        res["mean_mu"][i] = meanbranch(t -> LinearAlgebra.dot(model.μ, Ss[i](t)), t0, t1, x, w)
        res["delta_mu"][i] = LinearAlgebra.dot(model.μ, Ss[i](t0)) - LinearAlgebra.dot(model.μ, Ss[i](t1))

        ## net-diversification rate
        res["mean_netdiv"][i] = meanbranch(t -> LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t)), t0, t1, x, w)
        res["delta_netdiv"][i] = LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t0)) - LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t1))

        ## relative extinction rate (μ/λ)
        res["mean_relext"][i] = meanbranch(t -> LinearAlgebra.dot(model.μ ./ model.λ, Ss[i](t)), t0, t1, x, w)
        res["delta_relext"][i] = LinearAlgebra.dot(model.μ ./ model.λ, Ss[i](t0)) - LinearAlgebra.dot(model.μ ./ model.λ, Ss[i](t1))

        ## only if the model is an FBD model
        if hasproperty(model, :ψ)
            res["mean_psi"][i] = meanbranch(t -> LinearAlgebra.dot(model.ψ, Ss[i](t)), t0, t1, x, w)
            res["delta_psi"][i] = LinearAlgebra.dot(model.ψ, Ss[i](t0)) - LinearAlgebra.dot(model.ψ, Ss[i](t1))
        end
    end

    #edge = 1:n_branches
    edge = [branch.index for branch in branches]

    df = DataFrames.DataFrame(res)
    df[!, "node"] = get_node_indices(tree)
    df[!, "edge"] = edge

    tiplabs = tip_labels(tree)
    root_index = length(tiplabs)+1

    root_row = Dict(key => NaN for key in keys)
    root_row["node"] = root_index
    root_row["edge"] = 0

    df = push!(df, root_row)
    #push!(df, vcat(repeat([NaN], length(names)), root_index, 0))
    return(df)
end

#=
function tree_rates(tree::Root, model::T, Fs, Ss; n = 10) where {T <: ConstantModel}
    branches = get_branches(tree);
    n_branches = length(branches)

    rates = zeros(n_branches, 10)
    x, w = FastGaussQuadrature.gausslegendre(n)

    node_indices = Int64[]

    items = Vector{Float64}[]
    names = String[]
    
    #Threads.@threads for row in 1:size(data.edges)[1]
    for branch in branches
    #for branch in branches
        i = branch.index
        node_index = branch.outbounds.index
        push!(node_indices, node_index)
        t0, t1 = extrema(Fs[i].t)
        ## t0 is youngest, t1 is oldest

        ## posterior mean rate
        rates[i,1] = meanbranch(t -> LinearAlgebra.dot(model.λ, Ss[i](t)), t0, t1, x, w)
        rates[i,2] = meanbranch(t -> LinearAlgebra.dot(model.μ, Ss[i](t)), t0, t1, x, w)
        rates[i,3] = meanbranch(t -> LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t)), t0, t1, x, w)
        rates[i,4] = meanbranch(t -> LinearAlgebra.dot(model.μ ./ model.λ, Ss[i](t)), t0, t1, x, w)
        rates[i,5] = meanbranch(t -> LinearAlgebra.dot(model.ψ, Ss[i](t)), t0, t1, x, w)

        ## difference from oldest to youngest point on branch
        rates[i,6] = LinearAlgebra.dot(model.λ, Ss[i](t0)) - LinearAlgebra.dot(model.λ, Ss[i](t1))
        rates[i,7] = LinearAlgebra.dot(model.μ, Ss[i](t0)) - LinearAlgebra.dot(model.μ, Ss[i](t1))
        rates[i,8] = LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t0)) - LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t1))
        rates[i,9] = LinearAlgebra.dot(model.μ ./ model.λ, Ss[i](t0)) - LinearAlgebra.dot(model.μ ./ model.λ, Ss[i](t1))
        rates[i,10] = LinearAlgebra.dot(model.ψ, Ss[i](t0)) - LinearAlgebra.dot(model.ψ, Ss[i](t1))
    end
    #node = data.edges[:,2]
    #edge = 1:size(data.edges)[1]
   
    names = [
         "mean_lambda", "mean_mu", "mean_netdiv", "mean_relext", "mean_psi",
         "delta_lambda", "delta_mu", "delta_netdiv", "delta_relext", "delta_psi"
        ]
    df = DataFrames.DataFrame(rates, names)
    df[!, "node"] = node_indices
    df[!, "edge"] = collect(1:n_branches)
    root_index = tree.index
    push!(df, vcat(repeat([NaN], length(names)), root_index, 0))
    return(df)
end
=#

#=
function tree_rates(data::SSEdata, model::BDStimevarying, Fs, Ss; n = 10)
    rates = zeros(size(data.edges)[1], 8)
    x, w = FastGaussQuadrature.gausslegendre(n)
    
    Threads.@threads for i = 1:size(data.edges)[1]
    #for i = 1:size(data.edges)[1]
        t0, t1 = extrema(Fs[i].t)

        rates[i,1] = meanbranch(t -> LinearAlgebra.dot(model.λ(t), Ss[i](t)), t0, t1, x, w)
        rates[i,2] = meanbranch(t -> LinearAlgebra.dot(model.μ(t), Ss[i](t)), t0, t1, x, w)
        rates[i,3] = meanbranch(t -> LinearAlgebra.dot(model.λ(t) .- model.μ(t), Ss[i](t)), t0, t1, x, w)
        rates[i,4] = meanbranch(t -> LinearAlgebra.dot(model.μ(t) ./ model.λ(t), Ss[i](t)), t0, t1, x, w)

        ## difference from oldest to youngest point on branch
        ## t0 is youngest, t1 is oldest
        rates[i,5] = LinearAlgebra.dot(model.λ(t0), Ss[i](t0)) - LinearAlgebra.dot(model.λ(t1), Ss[i](t1))
        rates[i,6] = LinearAlgebra.dot(model.μ(t0), Ss[i](t0)) - LinearAlgebra.dot(model.μ(t1), Ss[i](t1))
        rates[i,7] = LinearAlgebra.dot(model.λ(t0) .- model.μ(t0), Ss[i](t0)) - LinearAlgebra.dot(model.λ(t1) .- model.μ(t1), Ss[i](t1))
        rates[i,8] = LinearAlgebra.dot(model.μ(t0)  ./ model.λ(t0), Ss[i](t0)) - LinearAlgebra.dot(model.μ(t1) ./ model.λ(t1), Ss[i](t1))
    end
    node = data.edges[:,2]
    edge = 1:size(data.edges)[1]
    names = ["mean_lambda", "mean_mu", "mean_netdiv", "mean_relext",
            "delta_lambda", "delta_mu", "delta_netdiv", "delta_relext"]
    df = DataFrames.DataFrame(rates, names)
    df[!, "node"] = node
    df[!, "edge"] = edge
    root_index = length(data.tiplab)+1
    push!(df, [NaN NaN NaN NaN NaN NaN NaN NaN root_index 0])
    return(df)
end
=#



function ancestral_state_probabilities(
        data::SSEdata, 
        post::Dict{Int64, OrdinaryDiffEq.ODESolution}, 
        pre::Dict{Int64, OrdinaryDiffEq.ODESolution}, 
    )
    Ss = Dict{Int64, Function}()
    for edge_idx in 1:(maximum(data.edges)-1)
        Ss[edge_idx] = t::Float64 -> pre[edge_idx](t) .* post[edge_idx](t)[:,2] ./ (sum(pre[edge_idx](t) .* post[edge_idx](t)[:,2]))
       #Ss[edge_idx] = t::Float64 -> Fs[edge_idx](t) .* Ds[edge_idx](t) ./ (sum(Fs[edge_idx](t) .* Ds[edge_idx](t)))
    end

    return (Ss)
end

function ancestral_state_probabilities(
        tree::Root, 
        post::Dict{Int64, OrdinaryDiffEq.ODESolution}, 
        pre::Dict{Int64, OrdinaryDiffEq.ODESolution}, 
    )
    branches = get_branches(tree)

    Ss = Dict{Int64,Function}()
    #for edge_idx in 1:(maximum(data.edges)-1)
    for branch in branches
        edge_idx = branch.index

        #Ss[edge_idx] = t::Float64 -> Fs[edge_idx](t) .* Ds[edge_idx](t) ./ (sum(Fs[edge_idx](t) .* Ds[edge_idx](t)))
        Ss[edge_idx] = t::Float64 -> pre[edge_idx](t) .* post[edge_idx](t)[:,2] ./ (sum(pre[edge_idx](t) .* post[edge_idx](t)[:,2]))
    end

    return (Ss)
end



## problem: this function is not type stable, or atleast S(t) is not 
function ancestral_state_probabilities(
        post::Dict{Int64, OrdinaryDiffEq.ODESolution}, 
        pre::Dict{Int64, OrdinaryDiffEq.ODESolution}, 
    )
    Ss = Dict{Int64,Function}()
    for edge_idx in collect(keys(post))
        Ss[edge_idx] = t::Float64 -> pre[edge_idx](t) .* post[edge_idx](t)[:,2] ./ (sum(pre[edge_idx](t) .* post[edge_idx](t)[:,2]))
    end

    return (Ss)
end

function ancestral_state_probability(D::Vector{Float64}, F::Vector{Float64}, t::Float64)
    S = D .* F ./ (sum(F .* D))
    return (S)
end
