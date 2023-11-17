export backwards_forwards_pass

function backwards_forwards_pass(model, data; alg = OrdinaryDiffEq.Tsit5())
    E = extinction_probability(model, data)
    Ds, sf = postorder(model, data, E; alg = alg)
    Fs = preorder(model, data, E, Ds; alg = alg)
    return(Ds, Fs)
end

