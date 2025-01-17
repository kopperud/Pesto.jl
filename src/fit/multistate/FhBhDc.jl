function FhBhDc_newmodel(x::Vector{T}; n = 6, sd = 0.587) where {T <: Real}
    #α, β = x[4:5]
    
    #rhat, ϵ, ψhat, α, β = x
    
    
    α = x[4]
    β = x[5]

    μhat = x[2] 
    λhat = μhat + maximum([x[1], α*5])
    ψhat = x[3]
   

    ## enforce that fossilization rate is bigger than fossilization shift rate?
    #ψmean = x[3] + β

    ## enforce that 
    # * lambda hat is bigger than mu
    # * speciation rate is bigger than shift rate in speciation

    
    #=
    dr = Distributions.LogNormal(log(rhat), sd)
    dψ = Distributions.LogNormal(log(ψhat), sd)
    
    rquantiles = Pesto.make_quantiles2(dr, n)
    ψ = Pesto.make_quantiles2(dψ, n)

    λ = rquantiles ./ (1 .- ϵ)
    μ = ϵ .* rquantiles ./ ( 1 .- ϵ )
    =#
    
    dλ = Distributions.LogNormal(log(λhat), sd)
    dψ = Distributions.LogNormal(log(ψhat), sd)

    λ = Pesto.make_quantiles2(dλ, n)
    μ = repeat([μhat], n)
    ψ = Pesto.make_quantiles2(dψ, n) 

    model = FhBhDcModel(λ, μ, ψ, α, β)

    return(model)
end

export fit_FhBhDc

function fit_FhBhDc(
    data::Root; 
    n = 6, 
    sd = 0.587, 
    n_attempts = 10,
    ## parameters    [    r,    μ,     ψ,    α,    β]
    lower          = [ 1e-05, 1e-5, 1e-04, 1e-8, 1e-8],
    upper          = [ 0.2, 0.35, 0.5, 0.3, 0.3],
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
        ps = [getpar(xi) for xi in x]
        #println("r: $(ps[1]) \t\t ϵ: $(ps[2]) \t ψ: $(ps[3]) \t α: $(ps[4]) \t β: $(ps[5])")
        println("λ: $(ps[1]+maximum([ps[2], 5*ps[4]])) \t\t μ: $(ps[2]) \t ψ: $(ps[3]) \t α: $(ps[4]) \t β: $(ps[5])")

        if any((x .- lower).^2 .< 1e-30)
            logl = -Inf
        elseif any((x .- upper) .^2 .< 1e-30)
            logl = -Inf
        else
            model = FhBhDc_newmodel(x; n = n, sd = sd)
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
        dβ = Distributions.LogNormal(log(0.01), 0.2)

        ## truncate the distribution
        tol = 1e-8
        dr = Distributions.Truncated(dr, lower[1] + tol, upper[1] - tol)
        dϵ = Distributions.Truncated(dϵ, lower[2] + tol, upper[2] - tol)
        dψ = Distributions.Truncated(dψ, lower[3] + tol, upper[3] - tol)
        dα = Distributions.Truncated(dα, lower[4] + tol, upper[4] - tol)
        dβ = Distributions.Truncated(dβ, lower[5] + tol, upper[5] - tol)
    else
        λml, μml = (0.3, 0.22)

        dλ = Distributions.LogNormal(log(0.5*λml), 0.5)
        dμ = Distributions.LogNormal(log(0.5*μml), 0.5)
        dψ = Distributions.LogNormal(log(0.05), 0.5)
        dα = Distributions.LogNormal(log(0.01), 0.2)
        dβ = Distributions.LogNormal(log(0.01), 0.2)

        tol = 1e-8
        dλ = Distributions.Truncated(dλ, lower[1] + tol, upper[1] - tol)
        dμ = Distributions.Truncated(dμ, lower[2] + tol, upper[2] - tol)
        dψ = Distributions.Truncated(dψ, lower[3] + tol, upper[3] - tol)
        dα = Distributions.Truncated(dα, lower[4] + tol, upper[4] - tol)
        dβ = Distributions.Truncated(dβ, lower[5] + tol, upper[5] - tol)
    end

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
            xinit[1] = rand(dλ)
            xinit[2] = rand(dμ)
            xinit[3] = rand(dψ)
            xinit[4] = rand(dα)
            xinit[5] = rand(dβ)
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
    model = FhBhDc_newmodel(x; n = n, sd = sd)
    return(optres, model, i-1)
end


