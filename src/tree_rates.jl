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
    I = quadrature(f, t0, t1, x, w)
    res = I / (t1 -  t0)
    return(res)
end

function tree_rates(data, model; n = 4)
    Ds, Fs = backwards_forwards_pass(model, data);
    Ss = ancestral_state_probabilities(data, Ds, Fs);

    tree_rates(data, model, Fs, Ss; n = n)
end

function tree_rates(data, model, Fs, Ss; n = 4)
    rates = zeros(size(data.edges)[1], 4)
    x, w = FastGaussQuadrature.gausslegendre(n)
    
    Threads.@threads for i = 1:size(data.edges)[1]
        t0, t1 = extrema(Fs[i].t)

        rates[i,1] = meanbranch(t -> LinearAlgebra.dot(model.λ, Ss[i](t)), t0, t1, x, w)
        rates[i,2] = meanbranch(t -> LinearAlgebra.dot(model.μ, Ss[i](t)), t0, t1, x, w)
        rates[i,3] = meanbranch(t -> LinearAlgebra.dot(model.λ .- model.μ, Ss[i](t)), t0, t1, x, w)
        rates[i,4] = meanbranch(t -> LinearAlgebra.dot(model.μ, Ss[i](t)) / LinearAlgebra.dot(model.λ, Ss[i](t)), t0, t1, x, w)
    end
    node_index = data.edges[:,2]
    edge_index = 1:size(data.edges)[1]
    names = ["mean_lambda", "mean_mu", "mean_netdiv", "mean_relext"]
    df = DataFrames.DataFrame(rates, names)
    df[!, "node_index"] = node_index
    df[!, "edge_index"] = edge_index
    root_index = length(data.tiplab)+1
    push!(df, [NaN NaN NaN NaN root_index 0])
    return(df)
end

function ancestral_state_probabilities(data, Ds, Fs)    
    Ss = Dict()
    for edge_idx in 1:(maximum(data.edges)-1)
       Ss[edge_idx] = t -> Fs[edge_idx](t) .* Ds[edge_idx](t) ./ (sum(Fs[edge_idx](t) .* Ds[edge_idx](t)))
    end

    return (Ss)
end

function calculate_tree_rates(data, model, Fs, Ss)

end

#= function calculate_tree_rates(data, model, Ds, Fs, Ss, nt = 100)
    average_branch_rates = Dict()
    for (rate, rate_name) in zip([model.λ, model.μ], ("λ", "μ"))
        d = zeros(size(data.edges)[1])

        Threads.@threads for i = 1:size(data.edges)[1]
            d[i] = average_branch_rate(Ss[i], Fs[i].t, rate, nt)
        end
        average_branch_rates[rate_name] = d
    end

    ancestors = make_ancestors(data)

    average_node_rates = Dict()
    for (rate, rate_name) in zip([model.λ, model.μ], ("λ", "μ"))
        m = zeros(maximum(data.edges))
        Threads.@threads for i in 1:maximum(data.edges)
            if i == length(data.tiplab)+1
                m[i] = NaN
            else
                edge_idx = ancestors[i]
                node_val = average_branch_rates[rate_name][edge_idx]
                m[i] = node_val
            end
        end
        average_node_rates[rate_name] = m
    end

    res = Dict("average_branch_rates" => average_branch_rates,
               "average_node_rates" => average_node_rates)

    return(res)
end =#
