export optimize_hyperparameters


function newmodel(x::Vector{T}; n = 6, H = 0.587) where {T <: Real}
    d, μmean, η = x
    λmean = d + μmean
    
    dλ = Distributions.LogNormal(log(λmean), H)
    dμ = Distributions.LogNormal(log(μmean), H)
    
    λquantiles = Pesto.make_quantiles2(dλ, n)
    µquantiles = Pesto.make_quantiles2(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)
    model = SSEconstant(λ, μ, η)

    return(model)
end

function getpar(x::Float64)
    return(x)
end

function getpar(x::ForwardDiff.Dual)
    return(x.value)
end

function optimize_hyperparameters(data::SSEdata; n = 6, H = 0.587)

    f(x::Vector{T}) where {T <: Real} = begin
        #d, μmean, η = x
        #λmean = d + μmean

        #println("λ: \t", getpar(λmean))
        #println("μ: \t", getpar(μmean))
        #println("η: \t", getpar(η))
        
        model = newmodel(x; n = n, H = H)

        logl = -logL_root(model, data)
        return(logl)
    end

    g!(G, x) = begin
        G[:] .= ForwardDiff.gradient(f, x)
    end


    inner_optimizer = Optim.GradientDescent()

    opts = Optim.Options(
            x_tol = 0.05, f_tol = 0.05, g_tol = 0.05, 
            show_trace = false, 
            iterations = 150, outer_iterations = 150)

    converged = false
    global i = 1
    while !converged && i <= 10
        ## random starting values
        xinit = rand(3) .* 0.1
        xinit[3] = xinit[3] * 0.1

        ## limits
        lower = [0.001, 0.001, 1e-08]
        upper = [2.0, 2.0, 0.8]

        global optres = Optim.optimize(f, g!, lower, upper, xinit, Optim.Fminbox(inner_optimizer), opts)

        converged = optres.x_converged || optres.f_converged || optres.g_converged
        i += 1

    end

    if !converged
        error("did not converge after $i iterations")
    end

    model = newmodel(optres.minimizer)
    return(optres, model)

end