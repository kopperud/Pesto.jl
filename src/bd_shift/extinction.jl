export extinction_probability

#notneg(u,p,t) = any(x->x<0,u)
isnegative(u,p,t) = any(x->x<0,u)
above_one(u,p,t) = any(x -> x >= 1, u)

function extinction_probability(model::Model, data::SSEdata) where {Model <: SSE}
    alg = OrdinaryDiffEq.Tsit5()
    K = number_of_states(model)
    #pE = (model.λ, model.μ, model.η, K)
    pE = (model, K)

    tree_height = maximum(data.node_depth)
    tspan = (0.0, tree_height)
    E0 = repeat([1.0 - data.sampling_probability], K)
    
    ode = extinction_prob(model)
    pr = OrdinaryDiffEq.ODEProblem{true}(ode, E0, tspan, pE);
    
    ## use low tolerance because we only solve E once, so we can afford it
    E = OrdinaryDiffEq.solve(pr, alg, abstol = 1e-10, reltol = 1e-10, isoutofdomain = above_one)#, isoutofdomain = ispositive)
    return(E)
end

function extinction_probability(model::Model, tree::Root)
    alg = OrdinaryDiffEq.Tsit5()
    K = number_of_states(model)
    pE = (model, K)

    tree_height = treeheight(tree)
    tspan = (0.0, tree_height)

    tip1 = find_one_extant_tip(tree)
    sampling_probability = tip1.sampling_probability

    E0 = repeat([1.0 - sampling_probability], K)
    
    ode = extinction_prob(model)
    pr = OrdinaryDiffEq.ODEProblem{true}(ode, E0, tspan, pE);
    
    ## use low tolerance because we only solve E once, so we can afford it
    E = OrdinaryDiffEq.solve(pr, alg, abstol = 1e-10, reltol = 1e-10, isoutofdomain = above_one)#, isoutofdomain = ispositive)
    
end

function extinction_probability(model::Model, sampling_probability::Float64, time::Float64)
    alg = OrdinaryDiffEq.Tsit5()
    K = number_of_states(model)
    p = (model, K)

    tspan = (0.0, time)

    E0 = repeat([1.0 - sampling_probability], K)
    
    ode = extinction_prob(model)
    pr = OrdinaryDiffEq.ODEProblem{true}(ode, E0, tspan, p);
    
    ## use low tolerance because we only solve E once, so we can afford it
    sol = OrdinaryDiffEq.solve(pr, alg, abstol = 1e-10, reltol = 1e-10, isoutofdomain = above_one, save_everystep = false)#, isoutofdomain = ispositive)
    E = sol.u[end] 
    return(E)
end

function extinction_probability(
        model::M, 
        sampling_probability::Float64, 
        time::Float64,
    ) where {M <: HomogeneousModel}

    alg = OrdinaryDiffEq.Tsit5()
    p = (model)
    tspan = (0.0, time)

    E0 = [1.0 - sampling_probability]
    
    ode = extinction_prob(model)
    pr = OrdinaryDiffEq.ODEProblem{true}(ode, E0, tspan, p);
    
    ## use low tolerance because we only solve E once, so we can afford it
    sol = OrdinaryDiffEq.solve(pr, alg, abstol = 1e-10, reltol = 1e-10, isoutofdomain = above_one, save_everystep = false)#, isoutofdomain = ispositive)
    E = sol.u[end][1]
    return(E)
end

