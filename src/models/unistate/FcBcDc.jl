export FcBcDcModel

struct FcBcDcModel{T <: Real} <: UniStateModel
    λ::T
    μ::T
    ψ::T
end

const FBDconstant = FcBcDcModel

function get_fossilization_rates(model::FcBcDcModel, t::Float64)
    return([model.ψ])
end

function get_speciation_rates(model::FcBcDcModel, t::Float64)
    return([model.λ])
end


function FcBcDc_ode(du, u, p, t)
    model, _ = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ

    du[1,1] = μ - (λ+μ+ψ)*u[1,1] + λ*u[1,1]*u[1,1] ## dE/dt
    du[1,2] = -(λ+μ+ψ)*u[1,2] + 2*λ*u[1,2]*u[1,1] ## dD/dt
end

function FcBcDc_forward_ode(du, u, p, t)
    model, _ = p
    λ = model.λ
    μ = model.μ
    ψ = model.ψ

    du[1,1] = - μ +(λ+μ+ψ)*u[1,1] - λ*u[1,1]*u[1,1] ## dE/dt
    du[1,2] = +(λ+μ+ψ)*u[1,2] - 2*λ*u[1,2]*u[1,1] ## dD/dt
end


function forward_prob(model::FcBcDcModel)
    return(FcBcDc_forward_ode)
end

function backward_prob(model::FcBcDcModel)
    return(FcBcDc_ode)
end


