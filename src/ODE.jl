## Probability that a lineage at time `t` is not represented in the reconstructed tree
## * This equation does not depend on the topology, so we solve for it first
function extinction_prob(dE, E, p, t)
    λ, μ, η, K = p
    dE[:] .= μ .- (λ .+ μ .+ η) .* E .+ λ .* E.^2 .+ (η/(K-1)) .* (sum(E) .- E) 
end

## Probability of of observing the branch at time `t`
## * We solve this equation in the postorder traversal
function backward_prob(dD, D, p, t)
    λ, μ, η, K, E = p

    Et = E(t)
    dD[:] .= - (λ .+ μ .+ η) .* D .+ 2 .* λ .* D .* Et .+ (η/(K-1)) .* (sum(D) .- D)
end

## This ODE is the previous one times minus one
## * We solve this equation in the preorder traversal, albeit with different starting values for each branch
function forward_prob(dF, F, p, t)
    λ, μ, η, K, E = p

    Et = E(t)
    dF[:] .= (-1) .* ( - (λ .+ μ .+ η) .* F .+ 2 .* λ .* F .* Et .+ (η/(K-1)) .* (sum(F) .- F))
end

## This is the ODE to solve for the numebr of rate shifts
function number_of_shifts!(dN, N, p, t)
    η, K, S, D = p

    Dt = D(t)
    St = S(t)
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

