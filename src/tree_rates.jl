export calculate_tree_rates

function calculate_tree_rates(data, model, Ds, Fs; verbose = false)
    if verbose
        println("Calculating state probabilities")
    end

    Ps = Dict()
    for edge_idx in 1:(maximum(data.edges)-1)
       Ps[edge_idx] = t -> Fs[edge_idx](t) .* Ds[edge_idx](t) ./ (sum(Fs[edge_idx](t) .* Ds[edge_idx](t)))
    end

    if verbose
        println("Calculating average branch rates")
    end

    average_branch_rates = Dict()
    for (rate, rate_name) in zip([model.λ, model.μ], ("λ", "μ"))
        d = Dict()
        for (key, P) in Ps
            times = Fs[key].t
            times1 = range(minimum(times), maximum(times), length = 100)
            d[key] = mean([sum(rate .* P) for P in Ps[key].(times1)])
        end
        average_branch_rates[rate_name] = d
    end

    if verbose
        println("Reordering for ape node indices")
    end

    average_node_rates = Dict()
    for (rate, rate_name) in zip([model.λ, model.μ], ("λ", "μ"))
        m = zeros(maximum(data.edges))
        for i in 1:maximum(data.edges)
            if i == length(data.tiplab)+1
                m[i] = NaN
            else
                edge_idx = findall(data.edges[:,2] .== i)[1]
                node_val = average_branch_rates[rate_name][edge_idx]
                m[i] = node_val
            end
        end
        average_node_rates[rate_name] = m
    end

    res = Dict("Ds" => Ds,
               "Fs" => Fs,
               "Ps" => Ps,
               "average_branch_rates" => average_branch_rates,
               "average_node_rates" => average_node_rates)
    return(res)
end

