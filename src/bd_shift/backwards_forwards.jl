export backwards_forwards_pass

function backwards_forwards_pass(
        model::Model, 
        data::SSEdata;
        condition = [:mrca, :survival],
    )

    #E = extinction_probability(model, data)
    Ds = postorder(model, data)
    Fs = preorder(model, data, Ds; condition = condition)
    return(Ds, Fs)
end

function backwards_forwards_pass(
        model::Model, 
        tree::Root;
        condition = [:mrca, :survival],
    )

    #E = extinction_probability(model, tree)
    Ds = postorder(model, tree)
    Fs = preorder(model, tree, Ds; condition = condition)

    return(Ds, Fs)
end

