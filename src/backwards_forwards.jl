export backwards_forwards_pass

function backwards_forwards_pass(model, data; verbose = false, alg = OrdinaryDiffEq.Tsit5())
    if verbose
        println("\t\t")
        println("Backwards-forwards method:")
    end

    E = extinction_probability(model, data)
    D_ends, Ds, sf = postorder(model, data, E; verbose = verbose, alg = alg)
    #logL = logL_root(model, data)
    #println("logL: \t", logL)

    Fs, F_ends = preorder(model, data, E, D_ends; verbose = verbose, alg = alg)

    return(Ds, Fs)
end

