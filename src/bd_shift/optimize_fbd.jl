function newmodel_fbd(x::Vector{T}; n = 6, sd = 0.587) where {T <: Real}
    α, β = x[4:5]

    μ = x[2]

    ## enforce that fossilization rate is bigger than fossilization shift rate?
    #ψmean = x[3] + β
    ψmean = x[3]

    ## enforce that 
    # * lambda hat is bigger than mu
    # * speciation rate is bigger than shift rate in speciation
    λmean = μ + x[1] + α  


    dλ = Distributions.LogNormal(log(λmean), sd)
    dψ = Distributions.LogNormal(log(ψmean), sd)
    
    λquantiles = Pesto.make_quantiles2(dλ, n)
    ψquantiles = Pesto.make_quantiles2(dψ, n)

    model = FBDSconstant(λmean, μ, ψmean, λquantiles, ψquantiles, α, β)

    return(model)
end

function optimize_hyperparameters2(
    data::Root; 
    n = 6, 
    sd = 0.587, 
    n_attempts = 10,
    lower = [1e-04, 1e-04, 1e-04, 1e-6, 1e-6],
    upper = [1.0, 1.0, 1.0, 0.3, 0.3],
    xinit = missing
    )

    ## create the logistic transform functions
    g, h = logistic(lower, upper, 0.5)

    f(x_tilde::Vector{T}) where {T <: Real} = begin
        #println([getpar(e) for e in x_tilde])
 
        x = g(x_tilde) ## backtransform to bounded realm
        #η = x[1]
        #μ = x[2]
        #λ = maximum([5*x[1], x[2]]) + x[3]
        #ps = [getpar(λ), getpar(μ), getpar(η)]
        #println("λ: $(ps[1]) \t\t μ: $(ps[2]) \t η: $(ps[3])")
        #println([getpar(e) for e in x])

        if any((x .- lower).^2 .< 1e-30)
            logl = -Inf
        elseif any((x .- upper) .^2 .< 1e-30)
            logl = -Inf
        else
            model = newmodel_fbd(x; n = n, sd = sd)
            logl = logL_root(model, data)
        end
        #println("logl: \t", getpar(logl))

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
    rml, μml = (0.1, 0.05)

    dμ = Distributions.LogNormal(log(0.5*μml), 0.5)
    dr = Distributions.LogNormal(log(0.5*rml), 0.5)
    dψ = Distributions.LogNormal(log(0.05), 0.5)

    dα = Distributions.LogNormal(log(0.01), 0.2)
    dβ = Distributions.LogNormal(log(0.01), 0.2)
    #dγ = Distributions.LogNormal(log(0.01), 0.5)

    ## truncate the distribution
    ϵ = 1e-8
    dμ = Distributions.Truncated(dμ, lower[1] + ϵ, upper[1] - ϵ)
    dr = Distributions.Truncated(dr, lower[2] + ϵ, upper[2] - ϵ)
    dψ = Distributions.Truncated(dψ, lower[3] + ϵ, upper[3] - ϵ)
    
    dα = Distributions.Truncated(dα, lower[4] + ϵ, upper[4] - ϵ)
    dβ = Distributions.Truncated(dβ, lower[5] + ϵ, upper[5] - ϵ)
    #dγ = Distributions.Truncated(dγ, lower[6] + ϵ, upper[6] - ϵ)

    inner_optimizer = Optim.Newton()

    opts = Optim.Options(
            #x_abstol = 0.05, f_abstol = 0.05, g_abstol = 0.05, 
            #x_tol = 0.05, f_tol = 0.05, g_tol = 0.05, 
            show_trace = false,
            iterations = 100, outer_iterations = 100)

    use_random_inits = ismissing(xinit)

    if use_random_inits
        xinit = zeros(5)
    end
    
    while !converged && i <= n_attempts

        if use_random_inits
            xinit[1] = rand(dμ)
            xinit[2] = rand(dr)
            xinit[3] = rand(dψ)

            xinit[4] = rand(dα)
            xinit[5] = rand(dβ)
            #xinit[6] = rand(dγ)
        end
            
        xinit_tilde = h(xinit)

        try
            global optres = Optim.optimize(f, g!, h!, xinit_tilde, inner_optimizer, opts)
            converged = optres.x_converged || optres.f_converged || optres.g_converged

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
    model = newmodel_fbd(x; n = n, sd = sd)
    return(optres, model, i-1)
end
