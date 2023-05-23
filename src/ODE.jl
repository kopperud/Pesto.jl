## Probability that a lineage at time `t` is not represented in the reconstructed tree
## * This equation does not depend on the topology, so we solve for it first
function extinction_prob_old(dE, E, p, t)
    λ, μ, η, i_not_js, K = p

    for i in 1:K
        i_not_j = i_not_js[i]
        dE[i] = μ[i] - (λ[i] + μ[i] + η) * E[i] + λ[i] * E[i]^2 + (η/(K-1)) * sum(E[j] for j in i_not_j) 
    end
end

## Vectorized ODE function
function extinction_prob(dE, E, p, t)
    λ, μ, η, K = p
    dE[:] .= μ .- (λ .+ μ .+ η) .* E .+ λ .* E.^2 .+ (η/(K-1)) .* (sum(E) .- E) 
end


## Probability of of observing the branch at time `t`
## * We solve this equation in the postorder traversal
function backward_prob_old(dD, D, p, t)
    λ, μ, η, i_not_js, K, E = p

    Et = E(t)
    for i in 1:K
        i_not_j = i_not_js[i]
        dD[i] = - (λ[i] + μ[i] + η) * D[i] + 2*λ[i]*D[i]*Et[i] + (η/(K-1)) * sum(D[j] for j in i_not_j)
    end
end
function backward_prob(dD, D, p, t)
    λ, μ, η, i_not_js, K, E = p

    Et = E(t)
    dD[:] .= - (λ .+ μ .+ η) .* D .+ 2 .* λ .* D .* Et .+ (η/(K-1)) .* (sum(D) .- D)
end

function backward_prob_outofplace(D, p, t)
    λ, μ, η, i_not_js, K, E = p

    Et = E(t)
    - (λ .+ μ .+ η) .* D .+ 2 .* λ .* D .* Et .+ (η/(K-1)) .* (sum(D) .- D)
end


## This ODE is the previous one times minus one
## * We solve this equation in the preorder traversal, albeit with different starting values for each branch
function forward_prob_old(dF, F, p, t)
    λ, μ, η, i_not_js, K, E = p

    Et = E(t)
    for i in 1:K
        i_not_j = i_not_js[i]
        dF[i] = (-1) * ( - (λ[i] + μ[i] + η) * F[i] + 2*λ[i]*F[i]*Et[i] + (η/(K-1)) * sum(F[j] for j in i_not_j))
    end
end

function forward_prob(dF, F, p, t)
    λ, μ, η, i_not_js, K, E = p

    Et = E(t)
    dF[:] .= (-1) .* ( - (λ .+ μ .+ η) .* F .+ 2 .* λ .* F .* Et .+ (η/(K-1)) .* (sum(F) .- F))
end


function number_of_shifts!(dN, N, p, t)
    η, K, S, D = p

    Dt = D(t)
    dN[:,:] .= - (ones(K) * ones(K)' .- LinearAlgebra.I(K)) .* (ones(K) * S(t)') .* (Dt * ones(K)') .* (η / (K-1)) ./ (ones(K) * Dt')
end