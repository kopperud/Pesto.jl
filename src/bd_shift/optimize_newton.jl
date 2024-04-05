export optimize_hyperparameters

export ConvergenceException

struct ConvergenceException <: Exception end

function optimize_hyperparameters(
    data::SSEdata; 
    n = 6, 
    sd = 0.587, 
    n_attempts = 20,
    #lower = [1e-08, 1e-04, 1e-04],
    upper = [0.4, 2.0, 1.0],
    xinit = missing
    )

    ## create the logistic transform functions
    g, h = logistic(upper, 0.5)

    f(x_tilde::Vector{T}) where {T <: Real} = begin
        x = g(x_tilde) ## backtransform to bounded realm

        η = x[1]
        μ = x[1] + x[2]
        λ = x[1] + x[2] + x[3]
        ps = [getpar(λ), getpar(μ), getpar(η)]
        #println("λ: $(ps[1]) \t\t μ: $(ps[2]) \t η: $(ps[3])")

        model = newmodel(x; n = n, sd = sd)

        
        logl = logL_root(model, data)

        # \t logl = $(getpar(logl))")

        #push!(history, [getpar(λmean), getpar(μmean), getpar(η), getpar(logl)])

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

    rml, μml = estimate_constant_netdiv_mu(data)

    dη = Distributions.LogNormal(log(0.01), 0.5)
    dμ = Distributions.LogNormal(log(0.5*μml), 0.5)
    dr = Distributions.LogNormal(log(0.5*rml), 0.5)

    ## truncate the distribution
    dη = Distributions.Truncated(dη, 0.0, upper[1])
    dμ = Distributions.Truncated(dμ, 0.0, upper[2])
    dr = Distributions.Truncated(dr, 0.0, upper[3])
    
    inner_optimizer = Optim.Newton()

    opts = Optim.Options(
            #x_abstol = 0.05, f_abstol = 0.05, g_abstol = 0.05, 
            #x_tol = 0.05, f_tol = 0.05, g_tol = 0.05, 
            show_trace = false,
            iterations = 150, outer_iterations = 150)
    
    while !converged && i <= n_attempts

        if ismissing(xinit)
            xinit = zeros(3)
            xinit[1] = rand(dη)
            xinit[2] = rand(dμ)
            xinit[3] = rand(dr)
        end
            

        println(xinit)
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
    model = newmodel(x; n = n, sd = sd)
    return(optres, model, i-1)
end