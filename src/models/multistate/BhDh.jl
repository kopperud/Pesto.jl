export BhDhModel, BDSconstant

struct BhDhModel{T1 <: Real, T2 <: Real} <: MultiStateModel
    λ::Vector{T1}
    μ::Vector{T1}
    η::T2
end

const BDSconstant = BhDhModel ## alias for old name
const SSEconstant = BhDhModel

function eltype(model::BhDhModel)
    return(typeof(model.η))
end

function get_speciation_rates(model::BhDhModel, t::Float64)
    return(model.λ)
end

function get_fossilization_rates(model::BhDhModel, time::Float64)
    error("can not get fossilization rate for a birth-death-shift model. either use an FBD model or don't include fossils in the tree.")
    #return(model.ψ)
end

function BhDh_extinction_ode(dE, E, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    η = model.η

    sumE = sum(E)

    LoopVectorization.@turbo warn_check_args=false for i in 1:K
        dE[i] = μ[i] - (λ[i] + μ[i] + η) * E[i] + λ[i] * E[i] * E[i] + (η/(K-1)) * (sumE - E[i]) 
    end
    nothing
end

function extinction_prob(model::BhDhModel)
    return(BhDh_extinction_ode)
end


## Probability of of observing the branch at time `t`
## * We solve this equation in the postorder traversal
function backward_ode(
        du::Matrix{T}, 
        u::Matrix{T}, 
        p, 
        t::Float64
    ) where {T <: Real}

    model, K = p
    E, D = eachcol(u)
    sumE = sum(E)
    sumD = sum(D)

    λ = model.λ
    μ = model.μ
    η = model.η

    r = η / (K-1)

    LoopVectorization.@turbo warn_check_args=false for i in axes(du, 1)
        du[i,1] = μ[i] -(λ[i]+μ[i]+η)*u[i,1] + λ[i]*u[i,1]*u[i,1] + r*(sumE-u[i,1]) 
        du[i,2] = -(λ[i]+μ[i]+η)*u[i,2] + 2*λ[i]*u[i,2]*u[i,1] + r*(sumD-u[i,2])
    end

    nothing
end

function backward_prob(model::BhDhModel)
    return(backward_ode)
end

## This ODE is the previous one times minus one
## * We solve this equation in the preorder traversal, albeit with different starting values for each branch
function BhDh_forward_ode(
        du::Matrix{T}, 
        u::Matrix{T}, 
        p, 
        t::Float64) where {T <: Real}

    model, K = p
    λ = model.λ
    μ = model.μ
    η = model.η

    sumE, sumF = sum(u, dims = 1)
    r = η / (K-1)

    #for i in axes(du, 1)
    LoopVectorization.@turbo warn_check_args=false for i in axes(du, 1)
        du[i,1] = - μ[i] + (λ[i] + μ[i] + η) * u[i,1] - λ[i] * u[i,1] * u[i,1] - r * (sumE - u[i,1]) 
        ## F
        du[i,2] = +(λ[i]+μ[i]+η)*u[i,2] - 2*λ[i]*u[i,2]*u[i,1] - r *(sumF - u[i,2])
    end
    nothing
end

function forward_prob(model::BhDhModel)
    return(BhDh_forward_ode)
end

## This is the ODE to solve for the number of rate shifts
function number_of_shifts!(dN, N, p, t)
    η, K, D, F = p

    Dt = D(t)[:,2]
    Ft = F(t)[:,2]
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

function shift_problem(model::BhDhModel)
    return(number_of_shifts!)
end

## this doesn't output a matrix but rather a scalar
function number_of_shifts_simple!(dN, N, p, t)
    η, K, D, F = p

    Dt = D(t)[:,2]
    Ft = F(t)[:,2]

    St = ancestral_state_probability(Dt, Ft, t)
    r = -(η/(K-1.0))

    ## [1] is the number of rate shifts, dN(t)/dt
    dN[1] = r * (sum(Dt .* sum(St ./ Dt)) - 1.0)
    nothing
end

function shift_problem_simple(model::BhDhModel)
    return(number_of_shifts_simple!)
end


function BhDh_no_shifts_prob(dlogX, logX, p, t)
    η, K, D = p

    r = η / (K-1.0)
    Dt = @view D(t)[:,2]
    Dsum = sum(Dt)

    # u[j] is the log probability that there were no shifts 
    # departing from state j from the beginning of the 
    # branch until time t.
    dlogX[:] = r .* (Dsum .- Dt) ./ Dt
end

function no_shifts_problem(model::BhDhModel)
    return(BhDh_no_shifts_prob)
end


