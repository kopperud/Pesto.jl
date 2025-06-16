export fit_BhDh

function BhDh_newmodel(x::Vector{T}; n = 6, sd = 0.587) where {T <: Real}
    η = x[1]
    μmean = x[2]
    λmean = maximum([5*x[1], x[2]]) + x[3]
            
    dλ = Distributions.LogNormal(log(λmean), sd)
    dμ = Distributions.LogNormal(log(μmean), sd)
    
    λquantiles = Pesto.make_quantiles2(dλ, n)
    µquantiles = Pesto.make_quantiles2(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)
    model = BDSconstant(λ, μ, η)

    return(model)
end


function fit_BhDh(
    data; 
    n = 6, 
    sd = 0.587, 
    n_attempts = 10,
    lower = [1e-08, 1e-04, 1e-04],
    upper = [0.3, 1.0, 1.0],
    xinit = missing,
    condition = [:mrca, :survival],
    ) 

    tl = tip_labels(data)
    ntips = length(tl)

    ## create the logistic transform functions
    g, h = logistic(lower, upper, 0.5)

    f(x_tilde::Vector{T}) where {T <: Real} = begin
        x = g(x_tilde) ## backtransform to bounded realm

        if any((x .- lower).^2 .< 1e-30)
            logl = -Inf
        elseif any((x .- upper) .^2 .< 1e-30)
            logl = -Inf
        else
            model = Pesto.BhDh_newmodel(x; n = n, sd = sd)
            logl = logL_root(model, data; condition = condition)
        end

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
    ϵ = 1e-8
    dη = Distributions.Truncated(dη, lower[1] + ϵ, upper[1] - ϵ)
    dμ = Distributions.Truncated(dμ, lower[2] + ϵ, upper[2] - ϵ)
    dr = Distributions.Truncated(dr, lower[3] + ϵ, upper[3] - ϵ)
    
    inner_optimizer = Optim.Newton()

    opts = Optim.Options(
            #x_abstol = 0.05, f_abstol = 0.05, g_abstol = 0.05, 
            #x_tol = 0.05, f_tol = 0.05, g_tol = 0.05, 
            show_trace = false,
            iterations = 100, outer_iterations = 100)

    use_random_inits = ismissing(xinit)

    if use_random_inits
        xinit = zeros(3)
    end
    
    while !converged && i <= n_attempts

        if use_random_inits
            xinit[1] = rand(dη)
            xinit[2] = rand(dμ)
            xinit[3] = rand(dr)
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
    model = BhDh_newmodel(x; n = n, sd = sd)
    return(optres, model, i-1)
end

const optimize_hyperparameters = fit_BhDh
