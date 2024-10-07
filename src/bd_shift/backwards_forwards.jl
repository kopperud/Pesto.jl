export backwards_forwards_pass

function backwards_forwards_pass(
        model::Model, 
        data::SSEdata; 
        alg = OrdinaryDiffEq.Tsit5())

    E = extinction_probability(model, data)
    Ds, sf = postorder(model, data, E; alg = alg)
    Fs = preorder(model, data, E, Ds; alg = alg)
    return(Ds, Fs)
end

function backwards_forwards_pass(model::Model, tree::Root) 

    E = extinction_probability(model, tree)
    Ds = postorder(model, tree, E)
    Fs = preorder(model, tree, E, Ds)

    return(Ds, Fs)
end

