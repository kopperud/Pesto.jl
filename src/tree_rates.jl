export calculate_tree_rates

function average_branch_rate(P, t, rate, nt)
    times = range(minimum(t), maximum(t), length = nt)

#    v = []
#    for time in times
#        append!(v, sum(rate .* P(time)))
#    end
    v = [sum(rate .* x) for x in P.(times)]
    res = StatsBase.mean(v)
    return(res)
end

function calculate_tree_rates(data, model, Ds, Fs; verbose = false, nt = 100)
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
        d = zeros(size(data.edges)[1])

        Threads.@threads for i = 1:size(data.edges)[1]
            d[i] = average_branch_rate(Ps[i], Fs[i].t, rate, nt)
        end
        average_branch_rates[rate_name] = d
    end

    if verbose
        println("Reordering for ape node indices")
    end

    ancestors = make_ancestors(data)

    average_node_rates = Dict()
    for (rate, rate_name) in zip([model.λ, model.μ], ("λ", "μ"))
        m = zeros(maximum(data.edges))
        Threads.@threads for i in 1:maximum(data.edges)
            if i == length(data.tiplab)+1
                m[i] = NaN
            else
    #            edge_idx = findall(data.edges[:,2] .== i)[1]
                edge_idx = ancestors[i]
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

    if verbose
        println("done ape reordering")
    end
    return(res)
end

#function average_node_rate(data, i, 


