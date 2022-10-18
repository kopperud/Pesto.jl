export backwards_forwards_pass

function backwards_forwards_pass(model, data; verbose = false)
    if verbose
        println("\t\t")
        println("Backwards-forwards method:")
    end

    D_ends, Ds, sf, E = postorder(model, data; verbose = verbose)
    #logL = logL_root(model, data)
    #println("logL: \t", logL)

    ASP, Fs = preorder(model, data, E, D_ends; verbose = verbose)
    if verbose
        display(ASP)
        println()
    end
    return (Ds, Fs)
end

