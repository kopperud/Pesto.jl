export lp, ψ, Econstant

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
    alg = DifferentialEquations.Tsit5()

    root_age = maximum(data.branching_times)
    tspan = (0.0, root_age)

    u0 = [E0, M0]
    p = [model.λ, model.μ]
    prob = DifferentialEquations.ODEProblem(dEdM, u0, tspan, p)
    EM = DifferentialEquations.solve(prob, alg, saveat = data.branching_times)
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

@doc raw"""
    loglikelihood(model, data)

From Louca and Pennell 2020 (Nature), eq. S28

```math
L = \frac{\rho^{n+1}} {\lambda (1 - E(t_1))^2} \times \prod_{i=1}^n \lambda \times \psi(0, t_i) \\
E(t) = 1 - \frac{\exp(\lambda - \mu)t}{\frac{1}{\rho} + \frac{\lambda}{\lambda -\mu} \Big ( \exp((\lambda - \mu)t) - 1 \Big)} \\
\psi(t) = \frac{e^{t(\lambda - \mu)}}{ [ 1 + \frac{\rho \lambda}{\lambda - \mu}(e^{t(\lambda - \mu)} - 1)]^{2}}
```

Logged:
```math
\log(L) = (n+1) \log(\rho) + \log(\psi(t_1)) - \log(\lambda) - 2 \log(1 - E(t_1)) + \sum_{i=1}^n \log(\lambda) + \log(\psi(t_i))
```

asd
"""
function lp(λ, μ, data::SSEdata)
#    λ = model.λ
#    μ = model.μ
    ρ = data.ρ

    ts = data.branching_times
    n = length(ts)

    logL = (n+1) * log(ρ) + log(ψ(ts[1], λ, μ, ρ))
    logL += - log(λ) - 2*log(1 - Econstant(ts[1], λ, μ, ρ))

    for i in 1:n
        logL += log(λ) + log(ψ(ts[i], λ, μ, ρ))
    end

    return(logL)
end

@doc raw"""
Equation S5 in Morlon et al. 2011 [PNAS]

```math
\psi(s, t) = e^{(\lambda - \mu)(t - s)} [ 1 + \frac{\frac{\lambda}{\lambda - \mu}(e^{t(\lambda - \mu)} - e^{s(\lambda-\mu)})}{\frac{1}{\rho} + \frac{\lambda}{\lambda - \mu} \times (e^{s(\lambda-\mu)}-1)}]^{-2}
```

We use this one, simplified where `s = 0`

```math
\psi(t) = \frac{e^{t(\lambda - \mu)}}{ [ 1 + \frac{\rho \lambda}{\lambda - \mu}(e^{t(\lambda - \mu)} - 1)]^{2}}
```
"""
function ψ(t, λ, μ, ρ)
    nom = exp(t * (λ - μ))
    denom = 1 + ((ρ * λ) /(λ - μ)) * (exp(t * (λ - μ)) - 1)
    res = nom / (denom*denom)

    return res
end

@doc raw"""
from Morlon et al. 2011 [PNAS], eq. S4

```math
E(t) = 1 - \frac{\exp(t(\lambda - \mu))}{\frac{1}{\rho} + \frac{\lambda}{\lambda -\mu} \Big ( \exp((\lambda - \mu)t) - 1 \Big)}
```
"""
function Econstant(t, λ, μ, ρ)
    nom = exp((λ - μ) * t)
    denom = (1 / ρ) + (λ / (λ - μ)) * (exp((λ - μ)*t) - 1)

    res = 1 - nom/denom
    return res
end

Turing.@model function birthdeath_constant(data)
    μ ~ Distributions.Exponential(0.1)
    d ~ Distributions.Exponential(0.1) # "net-diversification"

    λ = μ + d

    data ~ BDconstant(λ, μ)
end

function estimate_constant_bdp(data::SSEdata; iterations = 500)
    chain = Turing.sample(birthdeath_constant(data), Turing.NUTS(), iterations; progress=true)

    median_λ = StatsBase.median(chain[:d] + chain[:μ])
    median_μ = StatsBase.median(chain[:μ])
    return(chain, median_λ, median_μ)
end
