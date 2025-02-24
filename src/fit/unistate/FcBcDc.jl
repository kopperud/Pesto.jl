export fit_FcBcDc

function fit_FcBcDc(
        tree::Root; 
        xinit = [0.1, 0.05, 0.01], 
        lower = [0.000001, 0.000001, 0.000001], 
        upper = [20.0, 20.0, 10.0],
        condition = condition,
    )

    f(x_tilde) = begin
        x = exp.(x_tilde)
        model = FcBcDcModel(x[1], x[2], x[3])
        print("λ: ", getpar(x[1]), "\t μ: ", getpar(x[2]), "\t ψ: ", getpar(x[3]))
        lnl = logL_root(model, tree; condition = condition)
        println("\t logl: ", getpar(lnl))
        return(-lnl)
    end

    g!(G, x_tilde) = begin
        G[:] .= ForwardDiff.gradient(f, x_tilde)
    end

    ## updating the Hessian matrix
    h!(H, x_tilde) = begin
        H[:,:] = ForwardDiff.hessian(f, x_tilde)
    end

    inner_optimizer = Optim.Newton()
    xinit_tilde = log.(xinit)

    opts = Optim.Options(
            #x_abstol = 0.05, f_abstol = 0.05, g_abstol = 0.05, 
            #x_tol = 0.05, f_tol = 0.05, g_tol = 0.05, 
            show_trace = false,
            iterations = 100, outer_iterations = 100)

    optres = Optim.optimize(f, g!, h!, xinit_tilde, inner_optimizer, opts)
    λ, μ, ψ = exp.(optres.minimizer)
    m = FBDconstant(λ, μ, ψ)
    return(m)
end


