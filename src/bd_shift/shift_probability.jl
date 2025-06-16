export posterior_shift_prob
export prior_shift_prob
export posterior_prior_shift_odds



function no_shifts_prob_tv(dlogX, logX, p, t)
    η, K, D = p

    r = η(t) / (K-1.0)
    Dt = D(t)
    Dsum = sum(Dt)

    dlogX[:] = r .* (Dsum .- Dt) ./ Dt
end

#=
function no_shifts_problem(model::TimevaryingModel)
    return(no_shifts_prob_tv)
end
=#

isneg(u,p,t) = any(x->x>0,u)
#notneg(u,p,t) = any(x->x<0,u)
between0and1(u,p,t) = any(x -> (x < 0.0) || (x > 1.0), u)
#is_not_pos(u,p,t) = any(x -> x > 0, u)

function posterior_shift_prob_categories(model::Model, D, K, alg)
    t0 = D.t[1]
    t1 = D.t[end]
    tspan = (t1, t0)

    p = (model.η, K, D);
    logX0 = zeros(K);
    ode = no_shifts_problem(model)

    prob = OrdinaryDiffEq.ODEProblem(ode, logX0, tspan, p);

    sol = OrdinaryDiffEq.solve(
        prob, 
        alg, 
        #isoutofdomain = is_not_pos, ## log probs must 0 or smaller
        save_everystep = false);
    logX = sol.u[end]

    X = exp.(logX)

    return(X)
end

function posterior_shift_prob(model::Model, data::SSEdata)
    alg = OrdinaryDiffEq.Tsit5() 
    Ds, Fs = backwards_forwards_pass(model, data);

    n_edges = length(data.branch_lengths)

    K = number_of_states(model)
    X = zeros(n_edges)
    ode = no_shifts_problem(model)
    
    #Threads.@threads for edge_index in 1:n_edges
    for edge_index in 1:n_edges
        D = Ds[edge_index];
        F = Fs[edge_index];
        t1 = D.t[end]
    
        Xt0 = posterior_shift_prob_categories(model, D, K, alg)

        Dt = D(t1)[:,2]
        Ft = F(t1)#[:,2]
        St = ancestral_state_probability(Dt, Ft, t1)

        X[edge_index] = sum(Xt0 .* St)
    end
    prob_atleast_one_shift = 1.0 .- X
end

function posterior_shift_prob(model::MultiStateModel, tree::Root)
    alg = OrdinaryDiffEq.Tsit5() 
    Ds, Fs = backwards_forwards_pass(model, tree);

    #n_edges = length(tree.branch_lengths)
    n_branches = number_of_branches(tree)

    K = number_of_states(model)
    X = zeros(n_branches)
    ode = no_shifts_problem(model)
    
    Threads.@threads for edge_index in 1:n_branches
        D = Ds[edge_index];
        F = Fs[edge_index];
        t1 = D.t[end]
    
        Xt0 = posterior_shift_prob_categories(model, D, K, alg)
        St = ancestral_state_probability(D(t1), F(t1), t1)

        X[edge_index] = sum(Xt0 .* St)
    end
    prob_atleast_one_shift = 1.0 .- X
    
end



## https://en.wikipedia.org/wiki/Poisson_distribution
function poisson_pmf(model::BhDhModel, t0::Float64, t1::Float64, n::Int64)
    η = model.η
    time = t1 - t0
    r = η * time
    res = (r^n) * exp(-r) / factorial(n)
end

function poisson_zero(model::BhDhModel, t0::Float64, t1::Float64)
    η = model.η
    time = t1 - t0
    r = η * time
    res = exp(-r) 
end

#=
# https://gtribello.github.io/mathNET/resources/jim-chap22.pdf
function poisson_zero(model::TimevaryingModel, t0::Float64, t1::Float64)
    x, w = FastGaussQuadrature.gausslegendre(10)
    η_int = quadrature(model.η, t0, t1, x, w)
   
    res = exp(-η_int)
end
=#


function prior_shift_prob(model::Model, data::SSEdata)
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

function prior_shift_prob(model::Model, tree::Root)
    branches = get_branches(tree)
    n_branches = length(branches)

    prob_no_shift = zeros(n_branches)

    for branch in branches
        edge_index = branch.index

        bl = branch.time
        prob_no_shift[edge_index] = poisson_zero(model, 0.0, bl)
    end
    prob_shift = 1.0 .- prob_no_shift
    return(prob_shift)
end



# Shi, J. J., & Rabosky, D. L. (2015). Speciation dynamics during the global radiation of extant bats. Evolution, 69(6), 1528-1545.
function posterior_prior_shift_odds(model, data)
    prior_atleast_one_shift = prior_shift_prob(model, data)
    prior_no_shifts = 1.0 .- prior_atleast_one_shift

    posterior_atleast_one_shift = posterior_shift_prob(model, data)
    posterior_no_shifts = 1.0 .- posterior_atleast_one_shift

    odds = (posterior_atleast_one_shift ./ prior_atleast_one_shift) ./ (posterior_no_shifts ./ prior_no_shifts)
    return(odds)
end

