#export lp, psi, Econstant, estimate_constant_bdp, estimate_constant_fbdp

function printlambda(λ::Float64)
    println(λ)
end

function printlambda(λ::T) where {T <: ForwardDiff.Dual}
    println(λ.value)
end

#=
## the root
function postorder_async(
        model::M,
        root::Root,
    ) where {M <: UniStateModel}
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
        ) where {N <: BranchingEvent, M <: UniStateModel}

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
=#

#=
## along a branch
function postorder_async(
        model::M, 
        branch::Branch, 
        prob::OrdinaryDiffEq.ODEProblem,
        sampling_probability::Float64,
        time::Float64,
    ) where {M <: UniStateModel}
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
    ) where {M <: UniStateModel}

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
    ) where {M <: UniStateModel}
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
=#

#=
function logL_root(
        model::M, 
        tree::Root; 
        condition = [:survival, :mrca]
    ) where {M <: UniStateModel}
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
=#


