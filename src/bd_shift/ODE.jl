## Probability that a lineage at time `t` is not represented in the reconstructed tree
## * This equation does not depend on the topology, so we solve for it first
function extinction_ode(dE, E, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    η = model.η

    sumE = sum(E)

    dE[:] .= μ .- (λ .+ μ .+ η) .* E .+ λ .* E .* E .+ (η/(K-1)) .* (sumE .- E) 
end

function extinction_ode_tv(dE, E, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    η = model.η

    λ, μ, η, K = p
    dE[:] .= μ(t) .- (λ(t) .+ μ(t) .+ η(t)) .* E .+ λ(t) .* E .* E .+ (η(t)/(K-1)) .* (sum(E) .- E) 
end

function extinction_fossil_ode(dE, E, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    Q = model.Q
    #η = model.η
    K = number_of_states(model)

    #dE[:] .= μ .- (λ .+ μ .+ η .+ ψ) .* E .+ λ .* E .* E .+ (η/(K-1)) .* (sum(E) .- E) 
    dE[:] = μ .- (λ .+ μ .+ ψ) .* E .+ λ .* E .* E .+ Q * E 
end


function extinction_prob(model::BDSconstant)
    return(extinction_ode)
end
    
function extinction_prob(model::BDStimevarying)
    return(extinction_ode_tv)
end

function extinction_prob(model::FBDSconstant)
    return(extinction_fossil_ode)
end

## Probability of of observing the branch at time `t`
## * We solve this equation in the postorder traversal
function backward_ode(dD, D, p, t)
    model, K, E = p
    λ = model.λ
    μ = model.μ
    η = model.η

    Et = E(t)
    dD[:] .= - (λ .+ μ .+ η) .* D .+ 2 .* λ .* D .* Et .+ (η/(K-1)) .* (sum(D) .- D)
end


function backward_ode_tv(dD, D, p, t)
    λ, μ, η, K, E = p
    Et = E(t)
    dD[:] .= - (λ(t) .+ μ(t) .+ η(t)) .* D .+ 2 .* λ(t) .* D .* Et .+ (η(t)/(K-1)) .* (sum(D) .- D)
end

function backward_fossil_ode(dD, D, p, t)
    model, K, E = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    #η = model.η
    Q = model.Q

    Et = E(t)
    #dD[:] .= - (λ .+ μ .+ ψ .+ η) .* D .+ 2 .* λ .* D .* Et .+ (η/(K-1)) .* (sum(D) .- D)
    dD[:] .= - (λ .+ μ .+ ψ) .* D .+ 2 .* λ .* D .* Et .+ Q * D
end

function backward_fossil2_ode(dD, D, p, t)
    model, K, E = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    Q = model.Q

    Et = E(t)
    dD[:] = - (λ .+ μ .+ ψ) .* D .+ 2 .* λ .* D .* Et .+ Q * D 
end




function backward_prob(model::BDSconstant)
    return(backward_ode)
end

function backward_prob(model::BDStimevarying)
    return(backward_ode_tv)
end

function backward_prob(model::FBDSconstant)
    return(backward_fossil_ode)
end



## This ODE is the previous one times minus one
## * We solve this equation in the preorder traversal, albeit with different starting values for each branch
function forward_ode(dF, F, p, t)
    model, K, E = p
    λ = model.λ
    μ = model.μ
    η = model.η

    Et = E(t)
    dF[:] .= (-1) .* ( - (λ .+ μ .+ η) .* F .+ 2 .* λ .* F .* Et .+ (η/(K-1)) .* (sum(F) .- F))
end

function forward_fossil_ode(dF, F, p, t)
    model, K, E = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ
    #η = model.η
    Q = model.Q

    Et = E(t)
    #dF[:] .= (-1) .* ( - (λ .+ μ .+ ψ .+ η) .* F .+ 2 .* λ .* F .* Et .+ (η/(K-1)) .* (sum(F) .- F))
    dF[:] = (-1) .* ( - (λ .+ μ .+ ψ) .* F .+ 2 .* λ .* F .* Et .+ Q * F )
end



function forward_ode_tv(dF, F, p, t)
    λ, μ, η, K, E = p

    Et = E(t)
    dF[:] .= (-1) .* ( - (λ(t) .+ μ(t) .+ η(t)) .* F .+ 2 .* λ(t) .* F .* Et .+ (η(t)/(K-1)) .* (sum(F) .- F))
end

function forward_prob(model::BDSconstant)
    return(forward_ode)
end

function forward_prob(model::BDStimevarying)
    return(forward_ode_tv)
end

function forward_prob(model::FBDSconstant)
    return(forward_fossil_ode)
end

## this doesn't output a matrix but rather a scalar
function number_of_shifts_simple!(dN, N, p, t)
    η, K, D, F = p

    Dt = D(t)
    Ft = F(t)
    St = ancestral_state_probability(Dt, Ft, t)
    r = -(η/(K-1.0))

    ## [1] is the number of rate shifts, dN(t)/dt
    dN[1] = r * (sum(Dt .* sum(St ./ Dt)) - 1.0)
end

function number_of_shifts_simple_tv!(dN, N, p, t)
    η, K, D, F = p

    Dt = D(t)
    Ft = F(t)
    St = ancestral_state_probability(Dt, Ft, t)
    r = -(η(t)/(K-1.0))
 
    ## [1] is the number of rate shifts, dN(t)/dt
    dN[1] = r * (sum(Dt .* sum(St ./ Dt)) - 1.0)
end


## This is the ODE to solve for the number of rate shifts
function number_of_shifts!(dN, N, p, t)
    η, K, D, F = p

    Dt = D(t)
    Ft = F(t)
    St = ancestral_state_probability(Dt, Ft, t)
    r = -(η/(K-1.0))

    LoopVectorization.@turbo for i in 1:K, j in 1:K
        dN[i,j] = r * St[j] * Dt[i] / Dt[j]
    end
    ## assign diagonal zero afterwards, since LoopVectorization 
    ## does not know how to handle if statement
    LoopVectorization.@turbo for i in 1:K
        dN[i,i] = 0.0
    end
end

function number_of_shifts_tv!(dN, N, p, t)
    η, K, D, F = p

    Dt = D(t)
    Ft = F(t)
    St = ancestral_state_probability(Dt, Ft, t)
    r = -(η(t)/(K-1.0))

    LoopVectorization.@turbo for i in 1:K, j in 1:K
        dN[i,j] = r * St[j] * Dt[i] / Dt[j]
    end
    LoopVectorization.@turbo for i in 1:K
        dN[i,i] = 0.0
    end
end

function shift_problem(model::ConstantModel)
    return(number_of_shifts!)
end

function shift_problem(model::TimevaryingModel)
    return(number_of_shifts_tv!)
end

function shift_problem_simple(model::ConstantModel)
    return(number_of_shifts_simple!)
end

function shift_problem_simple(model::TimevaryingModel)
    return(number_of_shifts_simple_tv!)
end
