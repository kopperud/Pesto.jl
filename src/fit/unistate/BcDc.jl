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

phy = readtrees(Pesto.path("bears.tre"))
sampling_probability = 1.0
data = make_SSEdata(phy, sampling_probability)

lp(λ, μ, data)
```
"""
function lp(λ, μ, data::SSEdata)
    sampling_probability = data.sampling_probability

    ts = data.branching_times
    n = length(ts)

    logL = (n+1) * log(sampling_probability) + log(psi(ts[1], λ, μ, sampling_probability))
    logL += - log(λ) - 2*log(1 - Econstant(ts[1], λ, μ, sampling_probability))

    for i in 1:n
        logL += log(λ) + log(psi(ts[i], λ, μ, sampling_probability))
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
sampling_probability = 1.0
λ = 1.0
μ = 0.5
t = 0.1

psi(t, λ, μ, sampling_probability)
```
"""
function psi(t, λ, μ, sampling_probability)
    nom = exp(t * (λ - μ))
    denom = 1 + ((sampling_probability * λ) /(λ - μ)) * (exp(t * (λ - μ)) - 1)
    res = nom / (denom*denom)

    return res
end

@doc raw"""
from Morlon et al. 2011 [PNAS], eq. S4

```math
E(t) = 1 - \frac{\exp(t(\lambda - \mu))}{\frac{1}{\rho} + \frac{\lambda}{\lambda -\mu} \Big ( \exp((\lambda - \mu)t) - 1 \Big)}
```
"""
function Econstant(t, λ, μ, sampling_probability)
    nom = exp((λ - μ) * t)
    denom = (1 / sampling_probability) + (λ / (λ - μ)) * (exp((λ - μ)*t) - 1)

    res = 1 - nom/denom
    return res
end



@doc raw"""
    estimate_constant_bdp(data::SSEdata[; xinit = [0.11, 0.09], lower = [0.0001, 0.0001], upper = [20.0, 20.0]])

Estimates the speciation and extinction rate under the reconstructed birth-death process with time-homogeneous rates.

Example:

```julia
phy = readtree(Pesto.path("primates.tre"))
sampling_probability = 0.67
data = make_SSEdata(phy, sampling_probability)

λml, μml = estimate_constant_bdp(data)
```
"""
function estimate_constant_bdp(data::SSEdata; xinit = [0.11, 0.09], lower = [0.00000001, 0.00000001], upper = [20.0, 20.0])
    ## ML estimates of parameters
    f(x) = -lp(x[1], x[2], data) ## function to minimize
    g!(G, x) = begin
        G[:] .= ForwardDiff.gradient(f, x)
    end

    inner_optimizer = Optim.GradientDescent()
    optres = Optim.optimize(f, g!, lower, upper, xinit, Optim.Fminbox(inner_optimizer))

    λml, μml = optres.minimizer
    return(λml, μml)
end

function get_fossilization_rates(model::BcDcModel, t::Float64)
    error("birth-death model does not have a fossilization rate. either fit a fossilized birth-death model or use a tree without fossil samples.") 
end

function fit_BcDc(tree::Root; xinit = [0.11, 0.09], lower = [0.00000001, 0.00000001], upper = [20.0, 20.0])
    ## ML estimates of parameters
    f(x) = begin
        model = BcDcModel(x[1], x[2])
        return(-logL_root(model, tree))
    end

    g!(G, x) = begin
        G[:] .= ForwardDiff.gradient(f, x)
    end

    inner_optimizer = Optim.GradientDescent()
    optres = Optim.optimize(f, g!, lower, upper, xinit, Optim.Fminbox(inner_optimizer))

    λml, μml = optres.minimizer
    return(λml, μml)
end



export estimate_constant_netdiv_mu

function estimate_constant_netdiv_mu(data::SSEdata; xinit = [0.05, 0.1], lower = [0.00000001, 0.00000001], upper = [20.0, 20.0])
    ## ML estimates of parameters
    ## constrain λ > μ, i.e. r = λ - μ > 0

    ## x[1] is net-div
    ## x[2] is extinction rate
    ## x[2] + x[1] is speciation
    f(x) = begin  ## function to minimize
        return(-lp(x[2] + x[1], x[2], data))
    end
    g!(G, x) = begin
        G[:] .= ForwardDiff.gradient(f, x)
    end

    inner_optimizer = Optim.GradientDescent()
    optres = Optim.optimize(f, g!, lower, upper, xinit, Optim.Fminbox(inner_optimizer))

    r, μml = optres.minimizer
    #λml = μml + x2
    return(r, μml)
end

