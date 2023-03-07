export optimize_eta

@doc raw"""
    optimize_eta(λ, μ, data)

Finds the maximum-likelihood parameter value for η (the transition rate) under the birth-death-shift model with a finite state space, conditional on λ and μ.

Example:

```julia
using Pesto

phy = readtree(Pesto.path("primates.tre")) 
ρ = 0.67
data = make_SSEdata2(phy, ρ)

λ = [0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3] 
μ = [0.09, 0.09, 0.09, 0.19, 0.19, 0.19, 0.29, 0.29, 0.29] 

ηml = optimize_eta(λ, μ, data)

#model = SSEconstant(λ, μ, ηml)
```
11
"""
function optimize_eta(λ, μ, data; lower = 0.00001, upper = 10.0, xinit = 0.001)
    ## find the maximum-likelihood estimate of eta, the transition rate
    f(η) = -sselp(η[1], λ, μ, data)

    ## define the gradient function with respect to η
    g!(G, η) = begin
        G[1] = ForwardDiff.derivative(f, η[1])
    end

    #xinit = 1.0 / sum(data.branch_lengths)
    inner_optimizer = Optim.GradientDescent(linesearch=Optim.LineSearches.BackTracking(order=3))

    opts = Optim.Options(x_tol = 0.1, f_tol = 0.1, g_tol = 0.1, show_trace = false)
    result = Optim.optimize(f, g!, [lower], [upper], [xinit], Optim.Fminbox(inner_optimizer), opts)
    
    ηml = result.minimizer
    return(ηml)
end
