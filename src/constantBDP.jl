export lp, ψ, Econstant, estimate_constant_bdp

@doc raw"""
    lp(λ, μ, data)

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

Example:

```julia
λ = 1.0
μ = 0.5

phy = readtrees(Diversification.path("bears.tre"))
ρ = 1.0
data = make_SSEdata2(phy, ρ)

lp(λ, μ, data)
```
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

Example:

```julia
ρ = 1.0
λ = 1.0
μ = 0.5
t = 0.1

ψ(t, λ, μ, ρ)
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

@doc raw"""
    estimate_constant_bdp(data::SSEdata[; xinit = [0.11, 0.09], lower = [0.0001, 0.0001], upper = [20.0, 20.0]])

Estimates the speciation and extinction rate under the reconstructed birth-death process with time-homogeneous rates.

Example:

```julia
phy = readtrees(Diversification.path("bears.tre"))
ρ = 1.0
data = make_SSEdata2(phy, ρ)

λml, μml = estimate_constant_bdp(data)
```
"""
function estimate_constant_bdp(data::SSEdata; xinit = [0.11, 0.09], lower = [0.0001, 0.0001], upper = [20.0, 20.0])
    ρ = data.ρ

    ## ML estimates of parameters
    f(x) = -lp(x[1], x[2], data) ## function to minimize

    inner_optimizer = Optim.GradientDescent()
    optres = Optim.optimize(f, lower, upper, xinit, Fminbox(inner_optimizer))

    λml, μml = optres.minimizer
    return(λml, μml)
end
