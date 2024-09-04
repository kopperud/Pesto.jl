export extinction_probability

#notneg(u,p,t) = any(x->x<0,u)
isnegative(u,p,t) = any(x->x<0,u)
above_one(u,p,t) = any(x -> x >= 1, u)

function extinction_probability(model::Model, data::SSEdata) where {Model <: SSE}
    alg = OrdinaryDiffEq.Tsit5()
    K = number_of_states(model)
    pE = (model.λ, model.μ, model.η, K)

    tree_height = maximum(data.node_depth)
    tspan = (0.0, tree_height)
    E0 = repeat([1.0 - data.sampling_probability], K)
    
    ode = extinction_prob(model)
    pr = OrdinaryDiffEq.ODEProblem{true}(ode, E0, tspan, pE);
    
    ## use low tolerance because we only solve E once, so we can afford it
    E = OrdinaryDiffEq.solve(pr, alg, abstol = 1e-10, reltol = 1e-10, isoutofdomain = above_one)#, isoutofdomain = ispositive)
    return(E)
end
