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

function tree_rates(data::SSEdata, model::T, Fs, Ss; n = 10) where {T <: ConstantModel}
    rates = zeros(size(data.edges)[1], 8)
    x, w = FastGaussQuadrature.gausslegendre(n)
    
    Threads.@threads for i = 1:size(data.edges)[1]
    #for i = 1:size(data.edges)[1]
        t0, t1 = extrema(Fs[i].t)
        ## t0 is youngest, t1 is oldest

        ## posterior mean rate
        rates[i,1] = meanbranch(t -> LinearAlgebra.dot(model.λ, Ss[i](t)), t0, t1, x, w)
        rates[i,2] = meanbranch(t -> LinearAlgebra.dot(model.μ, Ss[i](t)), t0, t1, x, w)
        rates[i,3] = meanbranch(t -> LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t)), t0, t1, x, w)
        rates[i,4] = meanbranch(t -> LinearAlgebra.dot(model.μ ./ model.λ, Ss[i](t)), t0, t1, x, w)

        ## difference from oldest to youngest point on branch
        rates[i,5] = LinearAlgebra.dot(model.λ, Ss[i](t0)) - LinearAlgebra.dot(model.λ, Ss[i](t1))
        rates[i,6] = LinearAlgebra.dot(model.μ, Ss[i](t0)) - LinearAlgebra.dot(model.μ, Ss[i](t1))
        rates[i,7] = LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t0)) - LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t1))
        rates[i,8] = LinearAlgebra.dot(model.μ ./ model.λ, Ss[i](t0)) - LinearAlgebra.dot(model.μ ./ model.λ, Ss[i](t1))
    end
    node = data.edges[:,2]
    edge = 1:size(data.edges)[1]
    names = [
         "mean_lambda", "mean_mu", "mean_netdiv", "mean_relext",
         "delta_lambda", "delta_mu", "delta_netdiv", "delta_relext"
        ]
    df = DataFrames.DataFrame(rates, names)
    df[!, "node"] = node
    df[!, "edge"] = edge
    root_index = length(data.tiplab)+1
    push!(df, vcat(repeat([NaN], length(names)), root_index, 0))
    return(df)
end

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



function ancestral_state_probabilities(
        data::SSEdata, 
        post::Dict{Int64, OrdinaryDiffEq.ODESolution}, 
        pre::Dict{Int64, OrdinaryDiffEq.ODESolution}, 
    )
    Ss = Dict{Int64, Function}()
    for edge_idx in 1:(maximum(data.edges)-1)
        Ss[edge_idx] = t::Float64 -> pre[edge_idx](t)[:,2] .* post[edge_idx](t)[:,2] ./ (sum(pre[edge_idx](t)[:,2] .* post[edge_idx](t)[:,2]))
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
        Ss[edge_idx] = t::Float64 -> pre[edge_idx](t)[:,2] .* post[edge_idx](t)[:,2] ./ (sum(pre[edge_idx](t)[:,2] .* post[edge_idx](t)[:,2]))
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
        Ss[edge_idx] = t::Float64 -> pre[edge_idx](t)[:,2] .* post[edge_idx](t)[:,2] ./ (sum(pre[edge_idx](t)[:,2] .* post[edge_idx](t)[:,2]))
    end

    return (Ss)
end

function ancestral_state_probability(D::Vector{Float64}, F::Vector{Float64}, t::Float64)
    S = D .* F ./ (sum(F .* D))
    return (S)
end
