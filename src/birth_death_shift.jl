## Wrapper function
export birth_death_shift

function birth_death_shift(model, data; verbose = false)
    Ds, Fs = backwards_forwards_pass(model, data; verbose = verbose)
    res = calculate_tree_rates(data, model, Ds, Fs; verbose = verbose)

    average_node_rates = res["average_node_rates"]

    return(average_node_rates)
end
