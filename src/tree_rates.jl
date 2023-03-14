export calculate_tree_rates
export ancestral_state_probabilities

function average_branch_rate(S, t, rate, nt)
    times = range(minimum(t), maximum(t), length = nt)

#    v = zeros(nt)
#    for (i, time) in enumerate(times)
#        v[i] = sum(rate .* P(time))
#    end
    v = [sum(rate .* St) for St in S.(times)]
    res = Statistics.mean(v)
    return(res)
end

function ancestral_state_probabilities(data, model, Ds, Fs; verbose = false)
    if verbose
        println("Calculating state probabilities")
    end

    Ss = Dict()
    for edge_idx in 1:(maximum(data.edges)-1)
       Ss[edge_idx] = t -> Fs[edge_idx](t) .* Ds[edge_idx](t) ./ (sum(Fs[edge_idx](t) .* Ds[edge_idx](t)))
    end

    return (Ss)
end

function calculate_tree_rates(data, model, Ds, Fs, Ss, nt = 100)
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
end
