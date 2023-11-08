export backwards_forwards_pass

function backwards_forwards_pass(model, data; alg = OrdinaryDiffEq.Tsit5())
    E = extinction_probability(model, data)
    D_ends, Ds, sf = postorder(model, data, E; alg = alg)
    Fs, F_ends = preorder(model, data, E, D_ends; alg = alg)
    return(Ds, Fs)
end

