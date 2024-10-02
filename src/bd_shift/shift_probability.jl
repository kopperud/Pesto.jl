export posterior_shift_prob
export prior_shift_prob
export posterior_prior_shift_odds

function Qmatrix(model)
    K = length(model.λ)
    Q = zeros(K,K)
    Q .= model.η / (K-1)
    for i in 1:K
        Q[i,i] = - model.η
    end
    return(Q)
end

function Amatrix(model, E, Q, t)
    A = LinearAlgebra.diagm(- model.λ .- model.μ .+ 2 .* model.λ .* E(t)) .+ Q
    return(A) 
end


function P(model, t, Δt, D, E, Q)
    K = length(model.λ)
    Dt = D(t)
    A = Amatrix(model, E, Q, t)

    P1 = (LinearAlgebra.I(K) .- Δt .* A) .* (Dt * ones(K)')
    colSums = ones(K)*ones(K)' * P1

    res = P1 ./ colSums
    return(res)
end
#= 
function dumb_range(a::Float64, b::Float64, n::Int64)
    x = zeros(n)
    x[1] = a
    for i in 1:(n-1)
        x[]
    end
end =#

function posterior_shift_prob_difference_eq(model, data; n_knots = 20)
    E = extinction_probability(model, data);
    ## there is no point in factoring this out, because the rest of the function is much slower
    Ds, Fs = backwards_forwards_pass(model, data); 
    Ss = ancestral_state_probabilities(data, Ds, Fs);
    Q = Qmatrix(model)
    K = length(model.λ)
    n_edges = length(data.branch_lengths)

    prob_no_shift = zeros(n_edges, n_knots-1)
    times = [
        collect(range(Ds[i].t[1], Ds[i].t[end]; length = n_knots)) for i in 1:n_edges
    ]
    edge_indices = 1:n_edges
    Threads.@threads for edge_index in 1:n_edges
        t0 = Ds[edge_index].t[1]
        t1 = Ds[edge_index].t[end]
        span = t0 - t1
        Δt = span/(n_knots-1)
        S = Ss[edge_index]
        ts = times[edge_index]

        for i in 1:(n_knots-1)
            t = ts[i]
            St = S(t)
            Sm = ones(K) * transpose(St)
            P1 = P(model, t, Δt, Ds[edge_index], E, Q)
            P_ij = (1 .- LinearAlgebra.I(K)) .* P1 .* Sm
            P_no_shift = 1 .- sum(P_ij)
            prob_no_shift[edge_index,i] = P_no_shift
        end
    end

    shift_prob = 1 .- prod(prob_no_shift, dims = 2)[:,1]
    return(shift_prob)
end

#=
function no_shifts_prob(dlnX, lnX, p, t)
    D, F, η, K = p
    Dt = D(t)
    Ft = F(t)
    St = ancestral_state_probability(Dt, Ft, t)

    dlnX[1] = sum((η/(K-1)) .* St .* (sum(Dt) .- Dt) ./ Dt)
end

function no_shifts_prob_tv(dlnX, lnX, p, t)
    D, S, η, K = p
    Dt = D(t)
    dlnX[1] = sum((η(t)/(K-1)) .* S(t) .* (sum(Dt) .- Dt) ./ Dt)
end
=#

function no_shifts_prob(du, u, p, t)
    η, K, D = p

    r = η / (K-1.0)
    Dt = D(t)
    Dsum = sum(Dt)

    du[:] = r .* u .* (sum(Dt) .- Dt) ./ Dt
    #for i in 1:K
    #    du[i] = r * u[i] * (Dsum - Dt[i]) / Dt[i]
    #end
end

function no_shifts_prob_tv(du, u, p, t)
    η, K, D = p

    r = η(t) / (K-1.0)
    Dt = D(t)
    Dsum = sum(Dt)

    du[:] = r .* u .* (sum(Dt) .- Dt) ./ Dt
    #for i in 1:K
    #    du[i] = r * u[i] * (Dsum - Dt[i]) / Dt[i]
    #end
end



function no_shifts_problem(model::BDSconstant)
    return(no_shifts_prob)
end

function no_shifts_problem(model::BDStimevarying)
    return(no_shifts_prob_tv)
end

isneg(u,p,t) = any(x->x>0,u)

#=
function posterior_shift_prob(model::Model, data::SSEdata)
    alg = OrdinaryDiffEq.Tsit5()
    ## there is no point in factoring this out, because the rest of the function is much slower
    Ds, Fs = backwards_forwards_pass(model, data); 
    #Ss = ancestral_state_probabilities(data, Ds, Fs);

    n_edges = length(data.branch_lengths)
    K = number_of_states(model)
    lnX = zeros(n_edges)
    ode = no_shifts_problem(model)
    
    Threads.@threads for edge_index in 1:n_edges
        D = Ds[edge_index];
        F = Fs[edge_index];
        t0 = Ds[edge_index].t[1]
        t1 = Ds[edge_index].t[end]
        tspan = (t1, t0)

        p = (D, F, model.η, K);
        #u0 = zeros(1)
        u0 = Float64[0.0]
        prob = OrdinaryDiffEq.ODEProblem(ode, u0, tspan, p);

        sol = OrdinaryDiffEq.solve(
            prob, 
            alg, 
            isoutofdomain = isneg, ## log probabilities must be negative
            save_everystep = false);
        lnX[edge_index] = sol.u[end][1]
    end
    prob_atleast_one_shift = 1.0 .- exp.(lnX)
    return(prob_atleast_one_shift)
end
=#

#notneg(u,p,t) = any(x->x<0,u)
between0and1(u,p,t) = any(x -> (x < 0.0) || (x > 1.0), u)

function noshiftprob5_log(du, u, p, t)
    η, K, D = p

    r = η / (K - 1.0)
    Dt = D(t)

    du[:] .= r .* (sum(Dt) .- Dt) ./ Dt
end

#is_not_pos(u,p,t) = any(x -> x > 0, u)

function posterior_shift_prob_categories(model::BDSconstant, D, K, ode, alg)
    t0 = D.t[1]
    t1 = D.t[end]
    tspan = (t1, t0)

    p = (model.η, K, D);
    #u0 = ones(K);
    u0 = zeros(K);
    prob = OrdinaryDiffEq.ODEProblem(noshiftprob5_log, u0, tspan, p);

    sol = OrdinaryDiffEq.solve(
        prob, 
        alg, 
        #isoutofdomain = is_not_pos, ## log probs must 0 or smaller
        save_everystep = false);
    logX = sol.u[end]

    X = exp.(logX)
    # liitle hack here
    #X[X .> 1.0] .= 1.0

    return(X)
end

function posterior_shift_prob_categories2(model::Model, D, K, ode, alg)
    t0 = D.t[1]
    t1 = D.t[end]
    tspan = (t1, t0)

    p = (model.η, K, D);
    u0 = ones(K);
    prob = OrdinaryDiffEq.ODEProblem(ode, u0, tspan, p);

    sol = OrdinaryDiffEq.solve(
        prob, 
        alg, 
        #isoutofdomain = notneg, ## probabilities must be non-negative
        isoutofdomain = between0and1, ## probabilities must bt 0 and 1
        save_everystep = false);
    X = sol.u[end]

    # liitle hack here
    X[X .> 1.0] .= 1.0

    return(X)
end


function posterior_shift_prob(model::Model, data::SSEdata)
    alg = OrdinaryDiffEq.Tsit5() 
    Ds, Fs = backwards_forwards_pass(model, data);

    n_edges = length(data.branch_lengths)

    K = number_of_states(model)
    X = zeros(n_edges)
    ode = no_shifts_problem(model)
    
    Threads.@threads for edge_index in 1:n_edges
        D = Ds[edge_index];
        F = Fs[edge_index];
        t1 = D.t[end]
    
        Xt0 = posterior_shift_prob_categories(model, D, K, ode, alg)
        St = ancestral_state_probability(D(t1), F(t1), t1)

        #X[edge_index] = sum(sol.u[end] .* St)
        X[edge_index] = sum(Xt0 .* St)
    end
    prob_atleast_one_shift = 1.0 .- X
    
end

## https://en.wikipedia.org/wiki/Poisson_distribution
function poisson_pmf(model::BDSconstant, t0::Float64, t1::Float64, n::Int64)
    η = model.η
    time = t1 - t0
    r = η * time
    res = (r^n) * exp(-r) / factorial(n)
end

function poisson_zero(model::BDSconstant, t0::Float64, t1::Float64)
    η = model.η
    time = t1 - t0
    r = η * time
    res = exp(-r) 
end

# https://gtribello.github.io/mathNET/resources/jim-chap22.pdf
function poisson_zero(model::BDStimevarying, t0::Float64, t1::Float64)
    x, w = FastGaussQuadrature.gausslegendre(10)
    η_int = quadrature(model.η, t0, t1, x, w)
   
    res = exp(-η_int)
end


function prior_shift_prob(model, data)
    n_edges = length(data.branch_lengths)
    prob_no_shift = zeros(n_edges)

    for edge_index in 1:n_edges
        node_index = data.edges[edge_index,2]
        bl = data.branch_lengths[edge_index]
        t0 = data.node_depth[node_index]
        t1 = t0 + bl
        prob_no_shift[edge_index] = poisson_zero(model, t0, t1)
    end
    prob_shift = 1.0 .- prob_no_shift
    return(prob_shift)
end

# Shi, J. J., & Rabosky, D. L. (2015). Speciation dynamics during the global radiation of extant bats. Evolution, 69(6), 1528-1545.
function posterior_prior_shift_odds(model, data)
    prior_atleast_one_shift = Pesto.prior_shift_prob(model, data)
    prior_no_shifts = 1.0 .- prior_atleast_one_shift

    posterior_atleast_one_shift = posterior_shift_prob(model, data)
    posterior_no_shifts = 1.0 .- posterior_atleast_one_shift

    odds = (posterior_atleast_one_shift ./ prior_atleast_one_shift) ./ (posterior_no_shifts ./ prior_no_shifts)
    return(odds)
end

