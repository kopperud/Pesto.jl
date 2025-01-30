export calculate_equilibrium_frequencies

function ode_equilibrium_freq(
        dn::Vector{Float64}, 
        n::Vector{Float64}, 
        p, 
        t::Float64,
    )
    model, K = p;

    λ = model.λ
    μ = model.μ
    η = model.η

    r = 1 / (K-1)

    for j in axes(dn, 1)
        dn[j] = n[j] * (λ[j] - μ[j]) + η * (r*sum(n) - n[j]*(1 + r))
    end

    nothing
end

function calculate_equilibrium_frequencies(
        model::BhDhModel, 
        starting_frequencies::Vector{Float64},
        t::Float64,
    )

    n0 = starting_frequencies;
    tspan = (0.0, t)
    K = number_of_states(model)
    p = (model, K)

    prob = OrdinaryDiffEq.ODEProblem{true}(ode_equilibrium_freq, n0, tspan, p);
    sol = OrdinaryDiffEq.solve(prob)

    us = hcat(sol.u...)'
    
    freqs = us ./ sum(us, dims = 2)[:,1]

    return(sol.t, freqs)
end
