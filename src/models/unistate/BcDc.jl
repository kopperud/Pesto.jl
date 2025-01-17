export BcDcModel

struct BcDcModel{T <: Real} <: UniStateModel
    λ::T
    μ::T
end

const BDconstant = BcDcModel

## the ODE 
function BcDc_ode(du, u, p, t)
    model, _ = p
    λ = model.λ
    μ = model.μ

    du[1,1] = μ - (λ+μ)*u[1,1] + λ*u[1,1]*u[1,1] ## dE/dt
    du[1,2] = -(λ+μ)*u[1,2] + 2*λ*u[1,2]*u[1,1] ## dD/dt

    nothing
end

function get_speciation_rates(model::BcDcModel, t::Float64)
    return([model.λ]) 
end

function backward_prob(model::BcDcModel)
    return(BcDc_ode)
end

function number_of_states(model::BcDcModel)
    return (1)
end
