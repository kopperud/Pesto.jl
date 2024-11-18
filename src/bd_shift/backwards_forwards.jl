export backwards_forwards_pass

function backwards_forwards_pass(
        model::Model, 
        data::SSEdata,
    )

    #E = extinction_probability(model, data)
    Ds = postorder(model, data)
    Fs = preorder(model, data, Ds)
    return(Ds, Fs)
end

function backwards_forwards_pass(
        model::Model, 
        tree::Root,
    )

    #E = extinction_probability(model, tree)
    Ds = postorder(model, tree)
    Fs = preorder(model, tree, Ds)

    return(Ds, Fs)
end

