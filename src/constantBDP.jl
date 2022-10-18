function dEdM(du, u, p, t)
    λ, μ = p
    # dE/dt -- the extinction probability
    du[1] = μ - u[1] * (λ + μ) + λ * u[1]^2
    # dM/dt -- the number of expected lineages in the reconstructed tree
    du[2] = u[2] * λ * (u[1] - 1.0)
end

function Distributions.loglikelihood(model::BDconstant, data::SSEdata)
    E0 = 1.0 - data.ρ
    M0 = float(length(data.tiplab))
    alg = Tsit5()

    root_age = maximum(data.branching_times)
    tspan = (0.0, root_age)

    u0 = [E0, M0]
    p = [model.λ, model.μ]
    prob = ODEProblem(dEdM, u0, tspan, p)
    EM = solve(prob, alg, saveat = data.branching_times)
    E(t) = EM(t)[1]
    M(t) = EM(t)[2]

    dM(t) = EM(t, Val{1})[2] # The first derivative of M(t)
    n = length(data.branching_times)
    
    Mroot = M(root_age)
    if Mroot < 0
        return -Inf
    end
    logL = 2* log(M(root_age)) - (n+1) * log(M0)
    for i in 2:n
        x = -dM(data.branching_times[i])
        if x < 0
            return -Inf
        end
        logL += log(x)
    end
    return(logL)
end

using ForwardDiff
using Turing

@model function birthdeath_constant(data)
    μ ~ Exponential(0.1)
    d ~ Exponential(0.1) # "net-diversification"

    λ = μ + d

    data ~ BDconstant(λ, μ)
end

function estimate_constant_bdp(data::SSEdata; iterations = 500)
    chain = sample(birthdeath_constant(data), NUTS(), iterations; progress=true)

    median_λ = median(chain[:d] + chain[:μ])
    median_μ = median(chain[:μ])
    return(chain, median_λ, median_μ)
end
