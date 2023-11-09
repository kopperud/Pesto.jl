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

function posterior_shift_prob(model, data; n_knots = 20)
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

function no_shifts_prob(dlnX, lnX, p, t)
    D, S, η, K = p
    Dt = D(t)
    #dlnX[1] = sum((η/(K-1)) .* S(t) .* (1 .- D(t)) ./ D(t))
    dlnX[1] = sum((η/(K-1)) .* S(t) .* (sum(Dt) .- Dt) ./ Dt)
end


function posterior_shift_prob_ode(model, data)
    alg = OrdinaryDiffEq.Tsit5()
    ## there is no point in factoring this out, because the rest of the function is much slower
    Ds, Fs = backwards_forwards_pass(model, data); 
    Ss = ancestral_state_probabilities(data, Ds, Fs);

    n_edges = length(data.branch_lengths)
    K = length(model.λ)
    lnX = zeros(n_edges)
    
    Threads.@threads for edge_index in 1:n_edges
        S = Ss[edge_index]
        D = Ds[edge_index]
        t0 = Ds[edge_index].t[1]
        t1 = Ds[edge_index].t[end]
        tspan = (t1, t0)

        p = (D, S, model.η, K)
        u0 = zeros(1)
        prob = OrdinaryDiffEq.ODEProblem(no_shifts_prob, u0, tspan, p)
        sol = OrdinaryDiffEq.solve(prob, alg, save_everystep = false)
        lnX[edge_index] = sol[end][1]
    end
    prob_no_shift = 1 .- exp.(lnX)
    return(prob_no_shift)
end

## https://en.wikipedia.org/wiki/Poisson_distribution
function poisson_pmf(rate, time, n)
    r = rate * time
    res = (r^n) * exp(-r) / factorial(n)
end

function prior_shift_prob(model, data)
    n_edges = length(data.branch_lengths)
    prob_no_shift = zeros(n_edges)

    for edge_index in 1:n_edges
        time_span = data.branch_lengths[edge_index]
        prob_no_shift[edge_index] = poisson_pmf(model.η, time_span, 0)
    end
    prob_shift = 1.0 .- prob_no_shift
    return(prob_shift)
end

# Shi, J. J., & Rabosky, D. L. (2015). Speciation dynamics during the global radiation of extant bats. Evolution, 69(6), 1528-1545.
function posterior_prior_shift_odds(model, data; n_knots = 20)
    prior_atleast_one_shift = Pesto.prior_shift_prob(model, data)
    prior_no_shifts = 1.0 .- prior_atleast_one_shift
    #posterior_atleast_one_shift = posterior_shift_prob(model, data; n_knots = n_knots)
    posterior_atleast_one_shift = posterior_shift_prob_ode(model, data)
    posterior_no_shifts = 1.0 .- posterior_atleast_one_shift

    odds = (posterior_atleast_one_shift ./ prior_atleast_one_shift) ./ (posterior_no_shifts ./ prior_no_shifts)
    return(odds)
end

function bayes_factors(model, data)
    alg = OrdinaryDiffEq.tsit5()


end

