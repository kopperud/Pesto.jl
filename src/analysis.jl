## wrapper foranalysis

function pesto(data; n = 6, sd = 0.587)
    λml, μml = estimate_constant_bdp(data)

    dλ = Distributions.LogNormal(log(λml), sd)
    dμ = Distributions.LogNormal(log(μml), sd)

    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)
    η = optimize_eta(λ, µ, data)
    model = SSEconstant(λ, μ, η)

    Ds, Fs = backwards_forwards_pass(model, data);
    Ss = ancestral_state_probabilities(data, Ds, Fs);
    rates = calculate_tree_rates(data, model, Ds, Fs, Ss);
    nshift = compute_nshifts(model, data, Ds, Ss; ape_order = false)

    df = DataFrames.DataFrame(
        "edge_index" => 1:size(data.edges)[1],
        "nshift" => nshift,
        "mean_lambda" => rates["average_branch_rates"]["λ"],
        "mean_mu" => rates["average_branch_rates"]["μ"],
        "node_index" => data.edges[:,2],
    )
    df1 = DataFrames.DataFrame(
        "edge_index" => 0,
        "nshift" => 0,
        "mean_lambda" => 0,
        "mean_mu" => 0,
        "node_index" => length(data.tiplab)+1,
        )
    #df2 = hcat(df, df1)
    append!(df, df1)

    res = Dict(
        "treedata" => df,
        "model" => model,        
    )
    return(res)
end