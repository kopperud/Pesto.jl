export extinction_probability

function extinction_probability(model::Model, data::SSEdata) where {Model <: SSE}
    alg = OrdinaryDiffEq.Tsit5()
    K = number_of_states(model)
    pE = (model.λ, model.μ, model.η, K)

    tree_height = maximum(data.node_depth)
    tspan = (0.0, tree_height)
    E0 = repeat([1.0 - data.ρ], K)
    
    ode = extinction_prob(model)
    pr = OrdinaryDiffEq.ODEProblem{true}(ode, E0, tspan, pE);
    
    ## use low tolerance because we only solve E once, so we can afford it
    E = OrdinaryDiffEq.solve(pr, alg, abstol = 1e-10, reltol = 1e-10)
    return(E)
end
