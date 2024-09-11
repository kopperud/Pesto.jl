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

function posterior_variance(rate_categories::Vector{Float64}, St::Vector{Float64})#data::SSEdata, model::SSEconstant, Ss, t::Float64)
    m = LinearAlgebra.dot(rate_categories, St) 

    var = sum((rate_categories .- m) .^2 .* St)
    return(var)
end

function tree_rates(data::SSEdata, model::SSEconstant, Fs, Ss; n = 10)
    rates = zeros(size(data.edges)[1], 12)
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

        ## posterior variance
        rates[i, 5] = sqrt.(meanbranch(t -> posterior_variance(model.λ, Ss[i](t)), t0, t1, x, w))
        rates[i, 6] = sqrt.(meanbranch(t -> posterior_variance(model.μ, Ss[i](t)), t0, t1, x, w))
        rates[i, 7] = sqrt.(meanbranch(t -> posterior_variance(model.λ .- model.μ, Ss[i](t)), t0, t1, x, w))
        rates[i, 8] = sqrt.(meanbranch(t -> posterior_variance(model.μ ./ model.λ, Ss[i](t)), t0, t1, x, w))

        ## difference from oldest to youngest point on branch
        rates[i,9] = LinearAlgebra.dot(model.λ, Ss[i](t0)) - LinearAlgebra.dot(model.λ, Ss[i](t1))
        rates[i,10] = LinearAlgebra.dot(model.μ, Ss[i](t0)) - LinearAlgebra.dot(model.μ, Ss[i](t1))
        rates[i,11] = LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t0)) - LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t1))
        rates[i,12] = LinearAlgebra.dot(model.μ ./ model.λ, Ss[i](t0)) - LinearAlgebra.dot(model.μ ./ model.λ, Ss[i](t1))
    end
    node = data.edges[:,2]
    edge = 1:size(data.edges)[1]
    names = [
         "mean_lambda", "mean_mu", "mean_netdiv", "mean_relext",
         "std_lambda", "std_mu", "std_netdiv", "std_relext",
         "delta_lambda", "delta_mu", "delta_netdiv", "delta_relext"
        ]
    df = DataFrames.DataFrame(rates, names)
    df[!, "node"] = node
    df[!, "edge"] = edge
    root_index = length(data.tiplab)+1
    push!(df, vcat(repeat([NaN], length(names)), root_index, 0))
    return(df)
end

function tree_rates(data::SSEdata, model::SSEtimevarying, Fs, Ss; n = 10)
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



function ancestral_state_probabilities(data::SSEdata, Ds, Fs)    
    Ss = Dict()
    for edge_idx in 1:(maximum(data.edges)-1)
       Ss[edge_idx] = t::Float64 -> Fs[edge_idx](t) .* Ds[edge_idx](t) ./ (sum(Fs[edge_idx](t) .* Ds[edge_idx](t)))
    end

    return (Ss)
end

## problem: this function is not type stable, or atleast S(t) is not 
function ancestral_state_probabilities(Ds, Fs)    
    Ss = Dict()
    #for edge_idx in 1:(maximum(data.edges)-1)
    for edge_idx in collect(keys(Ds))
#       Ss[edge_idx] = t -> Fs[edge_idx](t) .* Ds[edge_idx](t) ./ (sum(Fs[edge_idx](t) .* Ds[edge_idx](t)))
       Ss[edge_idx] = t::Float64 -> Fs[edge_idx](t) .* Ds[edge_idx](t) ./ (sum(Fs[edge_idx](t) .* Ds[edge_idx](t)))
    end

    return (Ss)
end

function ancestral_state_probability(D::Vector{Float64}, F::Vector{Float64}, t::Float64)
    S = D .* F ./ (sum(F .* D))
    return (S)
end
