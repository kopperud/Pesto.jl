export lp, psi, Econstant, estimate_constant_bdp, estimate_constant_fbdp

function printlambda(λ::Float64)
    println(λ)
end

function printlambda(λ::T) where {T <: ForwardDiff.Dual}
    println(λ.value)
end


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

#    res = zeros(typeof(λ), n)
    #Threads.@threads for i in 1:n
    for i in 1:n
        logL += log(λ) + log(psi(ts[i], λ, μ, sampling_probability))
#        res[i] = log(λ) + log(psi(ts[i], λ, μ, sampling_probability))
    end
#   logL = sum(res)
#   logL += (n+1) * log(sampling_probability) + log(psi(ts[1], λ, μ, sampling_probability))
#   logL += - log(λ) - 2*log(1 - Econstant(ts[1], λ, μ, sampling_probability))

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

function estimate_constant_fbdp(tree::Root; xinit = [0.1, 0.05, 0.01], lower = [0.000001, 0.000001, 0.000001], upper = [20.0, 20.0, 10.0])

    f(x_tilde) = begin
        x = exp.(x_tilde)
        model = FBDconstant(x[1], x[2], x[3])
        print("λ: ", getpar(x[1]), "\t μ: ", getpar(x[2]), "\t ψ: ", getpar(x[3]))
        lnl = logL_root(model, tree)
        println("\t logl: ", getpar(lnl))
        return(-lnl)
    end

    g!(G, x_tilde) = begin
        G[:] .= ForwardDiff.gradient(f, x_tilde)
    end

    ## updating the Hessian matrix
    h!(H, x_tilde) = begin
        H[:,:] = ForwardDiff.hessian(f, x_tilde)
    end

    inner_optimizer = Optim.Newton()
    xinit_tilde = log.(xinit)

    opts = Optim.Options(
            #x_abstol = 0.05, f_abstol = 0.05, g_abstol = 0.05, 
            #x_tol = 0.05, f_tol = 0.05, g_tol = 0.05, 
            show_trace = false,
            iterations = 100, outer_iterations = 100)

    #optres = Optim.optimize(f, g!, lower, upper, xinit, Optim.Fminbox(inner_optimizer))
    optres = Optim.optimize(f, g!, h!, xinit_tilde, inner_optimizer, opts)
    λ, μ, ψ = exp.(optres.minimizer)
    m = FBDconstant(λ, μ, ψ)
    return(m)
end

## do it with tree iteration
function ode_constant_bdp(du, u, p, t)
    model = p
    λ = model.λ
    μ = model.μ

    du[1] = μ - (λ+μ)*u[1] + λ*u[1]*u[1] ## dE/dt
    du[2] = -(λ+μ)*u[2] + 2*λ*u[2]*u[1] ## dD/dt

    nothing
end

function ode_constant_fbdp(du, u, p, t)
    model = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ

    du[1] = μ - (λ+μ+ψ)*u[1] + λ*u[1]*u[1] ## dE/dt
    du[2] = -(λ+μ+ψ)*u[2] + 2*λ*u[2]*u[1] ## dD/dt
end

function backward_prob(model::FBDconstant)
    return(ode_constant_fbdp)
end

function backward_prob(model::BDconstant)
    return(ode_constant_bdp)
end

## the root
function postorder_async(
        model::M,
        root::Root,
    ) where {M <: HomogeneousModel}
    elt = eltype(model)
    K = number_of_states(model)
    #ode = backward_prob(model)
    ode = backward_prob(model)

    p = (model)
    tspan = (0.0, 1.0)
    u0 = ones(elt, 2)

    prob = OrdinaryDiffEq.ODEProblem{true,SciMLBase.FullSpecialize}(ode, u0, tspan, p)

    height = treeheight(root);

    ## find sampling probability
    ## assume it is equal in all the species
    leftmost_tip = find_one_extant_tip(root)
    sampling_probability = leftmost_tip.sampling_probability

    u, sf = postorder_async(model, root, prob, sampling_probability, height)
    return(u, sf)
end

## internal node
function postorder_async(
        model::M, 
        node::N, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64, 
        ) where {N <: BranchingEvent, M <: HomogeneousModel}

    branch_left, branch_right = node.children

    local u_left, u_right
    local sf_left, sf_right

    begin
        u_left, sf_left = postorder_async(model, branch_left, prob,sampling_probability,  time)
        u_right, sf_right = postorder_async(model, branch_right, prob,sampling_probability,  time)
    end

    D_left = u_left[2]
    D_right = u_right[2]
    E_left = u_left[1]

    D = D_left * D_right * model.λ
    c = sum(D)
    sf = sf_left + sf_right + log(c)
    D = D / c

    u = hcat(E_left, D)

    return(u, sf)
end


## along a branch
function postorder_async(
        model::M, 
        branch::Branch, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64,
    ) where {M <: HomogeneousModel}
    child_node = branch.outbounds

    t_old = time 
    t_young = time - branch.time

    u0, sf = postorder_async(model, child_node, prob, sampling_probability, t_young)

    tspan = (t_young, t_old)
    prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5(), isoutofdomain = notneg, save_everystep = false, reltol = 1e-3)

    u = sol.u[end]
    c = u[2]
    u[2] = u[2] ./ c

    if c > 0.0
        sf += log(c)
    else
        sf -= Inf
    end

    if !(sol.retcode == OrdinaryDiffEq.ReturnCode.Success)
        sf -= Inf
    end


    return(u, sf)
end


## for an extant tip
function postorder_async(
        model::M, 
        tip::ExtantTip, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64,
    ) where {M <: HomogeneousModel}

    @assert abs(time .- 0) < 0.001
    elt = eltype(model)
    K = number_of_states(model)

    E = ones(elt, K) .- tip.sampling_probability
    D = zeros(elt, K) .+ tip.sampling_probability
    sf = 0.0

    u = vcat(E, D)

    return(u, sf)
end

## for a fossil tip
function postorder_async(
        model::M,
        tip::FossilTip,
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64,
    ) where {M <: HomogeneousModel}
    elt = eltype(model)

    ## probability that sampled lineages immediately go extinct
    r = 0.0
    ψ = model.ψ

    Et = extinction_probability(model, sampling_probability, time)
    
    D = r * ψ + (1.0 - r) * ψ * Et
    sf = 0.0

    u = [Et, D]
    return(u, sf)
end

function logL_root(
        model::M, 
        tree::Root; 
        condition = [:survival, :mrca]
    ) where {M <: HomogeneousModel}
    u, sf = postorder_async(model, tree)

    E, D = u

    root_age = treeheight(tree)

    # we condition the likelihood by
    #
    # * that there was a speciation event at the MRCA
    # * that the two lineages subtending from the MRCA 
    #        must have survived until the present
    if :survival in condition
        nonextinct = (1.0 .- E).^2
        D = D ./ nonextinct
    end

    if :mrca in condition
        λroot = model.λ
        D = D ./ λroot
    end

    logL = log(D) + sum(sf)
    return(logL)
end


