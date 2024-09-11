## Probability that a lineage at time `t` is not represented in the reconstructed tree
## * This equation does not depend on the topology, so we solve for it first
function extinction_ode(dE, E, p, t)
    λ, μ, η, K = p
    dE[:] .= μ .- (λ .+ μ .+ η) .* E .+ λ .* E .* E .+ (η/(K-1)) .* (sum(E) .- E) 
end

function extinction_ode_tv(dE, E, p, t)
    λ, μ, η, K = p
    dE[:] .= μ(t) .- (λ(t) .+ μ(t) .+ η(t)) .* E .+ λ(t) .* E .* E .+ (η(t)/(K-1)) .* (sum(E) .- E) 
end

function extinction_prob(model::SSEconstant)
    return(extinction_ode)
end
    
function extinction_prob(model::SSEtimevarying)
    return(extinction_ode_tv)
end


## Probability of of observing the branch at time `t`
## * We solve this equation in the postorder traversal
function backward_ode(dD, D, p, t)
    λ, μ, η, K, E = p

    Et = E(t)
    dD[:] .= - (λ .+ μ .+ η) .* D .+ 2 .* λ .* D .* Et .+ (η/(K-1)) .* (sum(D) .- D)
end


function backward_ode_tv(dD, D, p, t)
    λ, μ, η, K, E = p
    Et = E(t)
    dD[:] .= - (λ(t) .+ μ(t) .+ η(t)) .* D .+ 2 .* λ(t) .* D .* Et .+ (η(t)/(K-1)) .* (sum(D) .- D)
end

function backward_prob(model::SSEconstant)
    return(backward_ode)
end

function backward_prob(model::SSEtimevarying)
    return(backward_ode_tv)
end


## This ODE is the previous one times minus one
## * We solve this equation in the preorder traversal, albeit with different starting values for each branch
function forward_ode(dF, F, p, t)
    λ, μ, η, K, E = p

    Et = E(t)
    dF[:] .= (-1) .* ( - (λ .+ μ .+ η) .* F .+ 2 .* λ .* F .* Et .+ (η/(K-1)) .* (sum(F) .- F))
end

function forward_ode_tv(dF, F, p, t)
    λ, μ, η, K, E = p

    Et = E(t)
    dF[:] .= (-1) .* ( - (λ(t) .+ μ(t) .+ η(t)) .* F .+ 2 .* λ(t) .* F .* Et .+ (η(t)/(K-1)) .* (sum(F) .- F))
end

function forward_prob(model::SSEconstant)
    return(forward_ode)
end

function forward_prob(model::SSEtimevarying)
    return(forward_ode_tv)
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

function shift_problem(model::SSEconstant)
    return(number_of_shifts!)
end

function shift_problem(model::SSEtimevarying)
    return(number_of_shifts_tv!)
end

function shift_problem_simple(model::SSEconstant)
    return(number_of_shifts_simple!)
end


function shift_problem_simple(model::SSEtimevarying)
    return(number_of_shifts_simple_tv!)
end
