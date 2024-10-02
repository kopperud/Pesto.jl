export optimize_gd


function newmodel(x::Vector{T}; n = 6, sd = 0.587) where {T <: Real}
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

function optimize_gd(
    data::SSEdata; 
    n = 6, 
    sd = 0.587, 
    n_attempts = 20,
    n_start_positions = 5)

    history = []



    f(x::Vector{T}) where {T <: Real} = begin
        η = x[1]
        μmean = x[1] + x[2]
        rmean = x[1] + x[2] + x[3]
        λmean = rmean + μmean
        
        model = newmodel(x; n = n, sd = sd)

        logl = logL_root(model, data)

        #println("λ: $(getpar(λmean)) \t\t μ: $(getpar(μmean)) \t η: $(getpar(η)) \t logl = $(getpar(logl))")

        push!(history, [getpar(λmean), getpar(μmean), getpar(η), getpar(logl)])

        return(-logl)
    end

    g!(G, x) = begin
        G[:] .= ForwardDiff.gradient(f, x)
        #println("gradient: \t", G)
    end


    inner_optimizer = Optim.GradientDescent()

    opts = Optim.Options(
            x_abstol = 0.05, f_abstol = 0.05, g_abstol = 0.05, 
            x_tol = 0.05, f_tol = 0.05, g_tol = 0.05, 
            show_trace = false, 
            iterations = 250, outer_iterations = 250)

    converged = false
    global i = 1

    rml, μml = estimate_constant_netdiv_mu(data)

    dr = Distributions.LogNormal(log(rml), 0.5)
    dμ = Distributions.LogNormal(log(μml), 0.5)
    dη = Distributions.LogNormal(log(0.01), 0.5)

    results = []
    global n_converged = 0

    while n_converged < n_start_positions && i <= n_attempts 
    
        #while !converged && i <= n_attempts # && n_converged < n_start_positions
        #while 

        ## random starting values
        #xinit = rand(3) .* 0.1
        #xinit[3] = xinit[3] * 0.1

        ## limits
        lower = [1e-08, 0.001, 0.001]
        upper = [0.8, 2.0, 1.0]

        xinit = zeros(3)
        for i in 1:100
            xinit[1] = rand(dη)
            xinit[2] = rand(dr)
            xinit[3] = rand(dμ)

            cond = all(xinit .> lower) & all(xinit .< upper)
            if cond
                continue
            end
        end
        if sum(xinit) == 0
            error("could not find initial values in boundaries")
        end

        global optres = Optim.optimize(f, g!, lower, upper, xinit, Optim.Fminbox(inner_optimizer), opts)

        converged = optres.x_converged || optres.f_converged || optres.g_converged

        if converged
            push!(results, optres)
            n_converged += 1
        end

        println("iter: $i, \t n_converged: $n_converged")

        i += 1
    end

    if !converged
        #error("did not converge $(n_converged) after $(i-1) iterations")
        println("did not converge $(n_converged) times after $(i-1) iterations")
    end
    
    return(results)

    model = newmodel(optres.minimizer; n = n, sd = sd)
    return(optres, model, i-1, history)

end
