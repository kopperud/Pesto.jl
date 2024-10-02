export optimize_eta

@doc raw"""
    optimize_eta(λ, μ, data)

Finds the maximum-likelihood parameter value for η (the transition rate) under the birth-death-shift model with a finite state space, conditional on λ and μ.

Example:

```julia
using Pesto

phy = readtree(Pesto.path("primates.tre")) 
sampling_probability = 0.635
data = SSEdata(phy, sampling_probability)

λ = [0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3] 
μ = [0.09, 0.09, 0.09, 0.19, 0.19, 0.19, 0.29, 0.29, 0.29] 

ηml = optimize_eta(λ, μ, data)

#model = BDSconstant(λ, μ, ηml)
```
11
"""
function optimize_eta(λ::Vector{Float64}, μ::Vector{Float64}, data; lower = -Inf, upper = Inf, xinit = -Inf)
    if !isfinite(xinit)
        xinit = minimum([
            0.1 / sum(data.branch_lengths),
            maximum(λ) / 4]
        )
    end

    if !isfinite(lower)
        #lower = 0.0001 * xinit
        lower = 1.0e-10
    end

    if !isfinite(upper)
        #upper = 100.0 * xinit
        upper = maximum(λ) / 2
    end

    ## find the maximum-likelihood estimate of eta, the transition rate
    f(η) = -sselp(η[1], λ, μ, data)

    ## define the gradient function with respect to η
    g!(G, η) = begin
        G[1] = ForwardDiff.derivative(f, η[1])
    end

#    inner_optimizer = Optim.GradientDescent(linesearch=Optim.LineSearches.BackTracking(order=3))
    inner_optimizer = Optim.GradientDescent()
#    inner_optimizer = Optim.Newton()    

    opts = Optim.Options(x_tol = 0.1, f_tol = 0.1, g_tol = 0.1, show_trace = false)
    result = Optim.optimize(f, g!, [lower], [upper], [xinit], Optim.Fminbox(inner_optimizer), opts)
    
    ηml = result.minimizer[1]
    return(ηml)
end

function optimize_eta(λ::Function, μ::Function, data; lower = -Inf, upper = Inf, xinit = -Inf)
    if !isfinite(xinit)
        xinit = minimum([
            0.1 / sum(data.branch_lengths),
            maximum(λ) / 4]
        )
    end

    if !isfinite(lower)
        #lower = 0.0001 * xinit
        lower = 1.0e-10
    end

    if !isfinite(upper)
        #upper = 100.0 * xinit
        upper = maximum(λ) / 2
    end

    ## find the maximum-likelihood estimate of eta, the transition rate
    f(η) = -sselp_tv(η[1], λ, μ, data)

    ## define the gradient function with respect to η
    g!(G, η) = begin
        G[1] = ForwardDiff.derivative(f, η[1])
    end

#    inner_optimizer = Optim.GradientDescent(linesearch=Optim.LineSearches.BackTracking(order=3))
    inner_optimizer = Optim.GradientDescent()
#    inner_optimizer = Optim.Newton()    

    opts = Optim.Options(x_tol = 0.1, f_tol = 0.1, g_tol = 0.1, show_trace = false)
    result = Optim.optimize(f, g!, [lower], [upper], [xinit], Optim.Fminbox(inner_optimizer), opts)
    
    ηml = result.minimizer[1]
    return(ηml)
end
