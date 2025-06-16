function FcBhDc_newmodel(x::Vector{T}; n = 6, sd = 0.587) where {T <: Real}
    α = x[4]

    μhat = x[2] 
    #λhat = μhat + x[1] + α
    λhat = x[1] + 5*α
    ψhat = x[3]
   
    ## enforce that 
    # * lambda hat is bigger than mu
    # * lambda is bigger than alpha

    dλ = Distributions.LogNormal(log(λhat), sd)
    λ = Pesto.make_quantiles2(dλ, n) 

    μ = repeat([μhat], n)
    ψ = repeat([ψhat], n)

    model = FcBhDcModel(λ, μ, ψ, α)

    return(model)
end

export fit_FcBhDc

function fit_FcBhDc(
    data::Root; 
    n = 6, 
    sd = 0.587, 
    n_attempts = 10,
    ## parameters    [    r,    μ,     ψ,    α]
    lower          = [ 1e-05, 1e-5, 1e-04, 1e-8],
    upper          = [ 0.2, 0.35, 0.5,     0.3],
    xinit = missing
    )

    ## create the logistic transform functions
    g, h = logistic(lower, upper, 0.5)

    f(x_tilde::Vector{T}) where {T <: Real} = begin
        #println([getpar(e) for e in x_tilde])
 
        x = g(x_tilde) ## backtransform to bounded realm
        ps = [getpar(xi) for xi in x]
        println("λ: $(ps[1]+maximum([ps[2], 5*ps[4]])) \t\t μ: $(ps[2]) \t ψ: $(ps[3]) \t α: $(ps[4])")

        if any((x .- lower).^2 .< 1e-30)
            logl = -Inf
        elseif any((x .- upper) .^2 .< 1e-30)
            logl = -Inf
        else
            model = FcBhDc_newmodel(x; n = n, sd = sd)
            logl = logL_root(model, data)
        end
        println("logl: \t", getpar(logl))

        return(-logl)
    end

    ## updating the gradient vector
    g!(G, x_tilde) = begin
        G[:] .= ForwardDiff.gradient(f, x_tilde)
    end

    ## updating the Hessian matrix
    h!(H, x_tilde) = begin
        H[:,:] = ForwardDiff.hessian(f, x_tilde)
    end

    converged = false
    global i = 1

    #rml, μml = estimate_constant_netdiv_mu(data)
    

    if false 
        rml, ϵml = (0.1, 0.6)

        dr = Distributions.LogNormal(log(0.5*rml), 0.5)
        dϵ = Distributions.LogNormal(log(0.5*ϵml), 0.5)
        dψ = Distributions.LogNormal(log(0.05), 0.5)
        dα = Distributions.LogNormal(log(0.01), 0.2)

        ## truncate the distribution
        tol = 1e-8
        dr = Distributions.Truncated(dr, lower[1] + tol, upper[1] - tol)
        dϵ = Distributions.Truncated(dϵ, lower[2] + tol, upper[2] - tol)
        dψ = Distributions.Truncated(dψ, lower[3] + tol, upper[3] - tol)
        dα = Distributions.Truncated(dα, lower[4] + tol, upper[4] - tol)
    else
        λml, μml = (0.3, 0.22)

        dλ = Distributions.LogNormal(log(0.5*λml), 0.5)
        dμ = Distributions.LogNormal(log(0.5*μml), 0.5)
        dψ = Distributions.LogNormal(log(0.05), 0.5)
        dα = Distributions.LogNormal(log(0.01), 0.2)

        tol = 1e-8
        dλ = Distributions.Truncated(dλ, lower[1] + tol, upper[1] - tol)
        dμ = Distributions.Truncated(dμ, lower[2] + tol, upper[2] - tol)
        dψ = Distributions.Truncated(dψ, lower[3] + tol, upper[3] - tol)
        dα = Distributions.Truncated(dα, lower[4] + tol, upper[4] - tol)
    end

    inner_optimizer = Optim.Newton()

    opts = Optim.Options(
            #x_abstol = 0.05, f_abstol = 0.05, g_abstol = 0.05, 
            #x_tol = 0.05, f_tol = 0.05, g_tol = 0.05, 
            show_trace = false,
            iterations = 100, outer_iterations = 100)

    use_random_inits = ismissing(xinit)

    if use_random_inits
        xinit = zeros(4)
    end
    
    while !converged && i <= n_attempts

        if use_random_inits
            xinit[1] = rand(dλ)
            xinit[2] = rand(dμ)
            xinit[3] = rand(dψ)
            xinit[4] = rand(dα)
        end
            
        xinit_tilde = h(xinit)

        try
            global optres = Optim.optimize(f, g!, h!, xinit_tilde, inner_optimizer, opts)
            converged = check_if_converged(optres)

        catch e
            if isa(e, AssertionError)
                println("assertion error with optimizing")
            else
                rethrow(e)
            end
        end


        i += 1
    end

    if !converged
        println("did not converge after $(i-1) iterations")
        throw(ConvergenceException())
        #println("did not converge $(n_converged) times after $(i-1) iterations")
    end
    
    x = g(optres.minimizer)
    model = FcBhDc_newmodel(x; n = n, sd = sd)
    return(optres, model, i-1)
end


