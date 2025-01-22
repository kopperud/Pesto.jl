## Probability that a lineage at time `t` is not represented in the reconstructed tree
## * This equation does not depend on the topology, so we solve for it first
function extinction_ode_matrix(dE, E, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    Q = model.Q

    #dE[:] .= μ .- (λ .+ μ .+ η) .* E .+ λ .* E .* E .+ (η/(K-1)) .* (sumE .- E) 
    LoopVectorization.@turbo warn_check_args=false for i in eachindex(E)
        x = μ[i] - (λ[i] + μ[i]) * E[i] + λ[i] * E[i] * E[i]

        for j in eachindex(E)
            ## this Q maybe it should be Q[j,i] 
            ## because we are going from young 
            ## to old, but since it is symmetric
            ## might be quicker to to Q[i,j]
            x += Q[i,j] * E[j] 
        end
        dE[i] = x
    end
    nothing

end

function extinction_ode_tv(dE, E, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    η = model.η

    dE[:] .= μ(t) .- (λ(t) .+ μ(t) .+ η(t)) .* E .+ λ(t) .* E .* E .+ (η(t)/(K-1)) .* (sum(E) .- E) 
end


function extinction_fossil_constant_ode(dE, E, p, t)
    model = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ

    dE[1] = μ - (λ+μ+ψ)*E[1] + λ*E[1]*E[1]
end

function extinction_prob(model::FBDconstant)
    return(extinction_fossil_constant_ode)
end



function backward_ode_matrix(du, u, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    Q = model.Q

    E, D = eachcol(u)
    dE, dD = eachcol(du)

    fastmv!(dE, Q, E)
    fastmv!(dD, Q, D)

    LoopVectorization.@turbo warn_check_args=false for i in axes(u, 1)
        du[i,1] += μ[i] -(λ[i]+μ[i])*u[i,1] + λ[i]*u[i,1]*u[i,1] 
        du[i,2] += -(λ[i]+μ[i])*u[i,2] + 2*λ[i]*u[i,2]*u[i,1]
    end

    nothing
end

function forward_ode_matrix(du, u, p, t)
    backward_ode_matrix(du, u, p, t)

    du[:,:] = (-1) .* du[:,:]
    nothing
end

## computes the product A*x
## and stores into y
function fastmv!(y, A, x)
    LoopVectorization.@turbo warn_check_args=false for i in axes(A, 1)
        yi = zero(Base.eltype(x))
        for j in axes(A, 2)
            yi += A[i,j] * x[j]
        end
        y[i] = yi        
    end
    nothing
end



#=
function backward_ode_tv(dD, D, p, t)
    λ, μ, η, K, E = p
    Et = E(t)
    dD[:] .= - (λ(t) .+ μ(t) .+ η(t)) .* D .+ 2 .* λ(t) .* D .* Et .+ (η(t)/(K-1)) .* (sum(D) .- D)
end
=#

#=
function backward_prob(model::BDSconstantQ)
    return(backward_ode_matrix)
end
=#

#=
function backward_prob(model::BDStimevarying)
    return(backward_ode_tv)
end
=#

#=
function forward_ode_tv(dF, F, p, t)
    λ, μ, η, K, E = p

    Et = E(t)
    dF[:] .= (-1) .* ( - (λ(t) .+ μ(t) .+ η(t)) .* F .+ 2 .* λ(t) .* F .* Et .+ (η(t)/(K-1)) .* (sum(F) .- F))
end
=#

#=
function forward_prob(model::BDStimevarying)
    return(forward_ode_tv)
end
=#


#=
function forward_prob(model::BDSconstantQ)
    return(forward_ode_matrix)
end
=#

#=
function number_of_shifts_simple_tv!(dN, N, p, t)
    η, K, D, F = p

    Dt = D(t)[:,2]
    Ft = F(t)[:,2]

    St = ancestral_state_probability(Dt, Ft, t)
    r = -(η(t)/(K-1.0))
 
    ## [1] is the number of rate shifts, dN(t)/dt
    dN[1] = r * (sum(Dt .* sum(St ./ Dt)) - 1.0)
end


function number_of_shifts_tv!(dN, N, p, t)
    η, K, D, F = p

    Dt = D(t)[:,2]
    Ft = F(t)[:,2]
    St = ancestral_state_probability(Dt, Ft, t)
    r = -(η(t)/(K-1.0))

    LoopVectorization.@turbo for i in 1:K, j in 1:K
        dN[i,j] = r * St[j] * Dt[i] / Dt[j]
    end
    LoopVectorization.@turbo for i in 1:K
        dN[i,i] = 0.0
    end
end
=#

#=
function shift_problem(model::TimevaryingModel)
    return(number_of_shifts_tv!)
end
=#


#=
function shift_problem_simple(model::TimevaryingModel)
    return(number_of_shifts_simple_tv!)
end
=#
