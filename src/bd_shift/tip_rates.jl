export tip_rates

function tip_rates(model::SSEconstant, data::SSEdata, Ds, Fs)
    Ss = ancestral_state_probabilities(data, Ds, Fs)
    ntips = length(data.tiplab)
    ancestors = make_ancestors(data)

    res = zeros(ntips, 4)
   
    rates = (
        model.λ,
        model.µ,
        model.λ .- model.µ,
        model.µ ./ model.λ
    )
   
    for i in 1:size(data.tiplab)[1]
        edge_index = ancestors[i]
        for (j, rate) in enumerate(rates)
            t = Ds[edge_index].t[1]
            S = Ss[edge_index](t)

            res[i,j] = rate' * S
        end        
    end

    names = ["lambda", "mu", "netdiv", "relext"]
    df = DataFrames.DataFrame(res, names)
    df[!,"species"] = data.tiplab

    return(df)
end

function tip_rates(model::SSEconstant, data::SSEdata)
    Ds, Fs = backwards_forwards_pass(model, data)

    tip_rates(model, data, Ds, Fs)
end