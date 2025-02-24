export branch_variance

export sp_mean
export sp_variance

function sp_mean(model, Ss, branch_index, t)
    m = sum(model.λ .* Ss[branch_index](t))

    return(m)
end

function sp_variance(model, Ss, branch_index, t)
    # 1st moment
    m1 = sum(model.λ .* Ss[branch_index](t))

    # 2nd moment
    m2 = sum(model.λ .^ 2 .* Ss[branch_index](t))

    v = m2 - m1^2

    return(v)
end

#=
export Amatrix, Pmatrix

function Amatrix(model, E, K, t)
    Q = -LinearAlgebra.I(K) .* model.η .+ (1 .- LinearAlgebra.I(K)) .* (model.η/(K-1))
    A = LinearAlgebra.diagm(- (model.λ .+ model.μ) .+ 2 .* model.λ .* E(t)) .+ Q
    return(A)
end

function Pmatrix(model, D, E, t, Δt)
    K = length(model.λ)
    A = Amatrix(model, E, K, t)
    Dt = D(t+Δt)[:,2]

    P_unnorm = (LinearAlgebra.I(K) .- Δt .* A) .* (Dt * ones(K)')
    csum = ones(K) * sum(P_unnorm, dims = 1) ## column sum
    P = P_unnorm ./ csum
    return(P)
end
=#

function branch_variance(model, data, c)
    Ds, Fs = backwards_forwards_pass(model, data)

    vs = zeros(Float64, number_of_branches(data))
    ProgressMeter.@showprogress for edge_index in eachindex(Ds)
        D = Ds[edge_index]
        F = Fs[edge_index]

        v = branch_variance(model, data, D, F, c)
        vs[edge_index] = v
    end
    return(vs)
end

function branch_variance(
        model::BhDhModel,
        data, 
        branch_index::Int64,
        c::Int64)
    Ds, Fs = backwards_forwards_pass(model, data)
    F = Fs[branch_index]
    D = Ds[branch_index]

    v = branch_variance(model, data, D, F, c)
    return(v)
end


function branch_variance(
        model::BhDhModel, 
        data, 
        D::OrdinaryDiffEq.ODESolution, 
        F::OrdinaryDiffEq.ODESolution, 
        c::Int64)

    domain = D.t
    tmin = minimum(domain)
    tmax = maximum(domain)

    times = collect(range(tmax, tmin; length = c))
    
    var_barlambda = 0.0
    times = collect(range(D.t[end], D.t[1]; length = c))

    for a in 1:c
        t = times[a]
        
        Dt = D(t)[:,2]
        Ft = F(t)[:,2]
        St = Dt .* Ft
        St = St / sum(St)

        first_moment = sum(model.λ .* St)
        second_moment = sum(model.λ .* model.λ .* St)

        va = second_moment - first_moment^2
        var_barlambda += (1/c^2) * va

        for b in 2:c
            if b > a
                t_b = times[b]
                var_barlambda += 2*(1/c^2) * covariance(t, t_b, model, D, F)
            end
        end
    end
    return(var_barlambda)
end



function covariance(
        t1::Float64, 
        t2::Float64,
        model::BhDhModel, 
        D::OrdinaryDiffEq.ODESolution, 
        F::OrdinaryDiffEq.ODESolution, 
    )
    D1 = D(t1)[:,2]
    D2 = D(t2)[:,2]

    F1 = F(t1)[:,2]
    F2 = F(t2)[:,2]

    S1 = D1 .* F1
    S1 = S1 / sum(S1)

    S2 = D2 .* F2
    S2 = S2 / sum(S2)

    P = Pesto.Pmatrix_precise(model, D, t1, t2)

    E1 = sum(model.λ .* S1) ## expectation at time t1
    E2 = sum(model.λ .* S2) ## expectation at time t2

    ## covariance
    co = 0.0
    for i in 1:4
        for j in 1:4
            co += S1[j] * P[i,j] * (model.λ[j] - E1) * (model.λ[i] - E2)
        end
    end
    return(co)
end


function Pmatrix_precise(model, D, t1, t2)
    K = number_of_states(model)
    P = zeros(K,K)
    tspan = (t1, t2)
    p = (model, K)
    Et = D(t1)[:,1]

    D2 = D(t2)[:,2]

    for j in 1:K
        u0 = zeros(K,2)
        u0[:,1] = Et
        u0[j,2] = 1.0
        ode = BhDh_forward_ode

        prob = OrdinaryDiffEq.ODEProblem{true}(ode, u0, tspan, p)
        sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5(), isoutofdomain = notneg, save_everystep = false, reltol = 1e-12)
        F2 = sol.u[end][:,2] 
        
        P[:,j] = F2 .* D2
        #P[:,j] = F2 
    end

    csum = ones(K) * sum(P, dims = 1)
    P = P ./ csum
    

    #rsum = sum(P, dims = 2) * ones(K)'
    #P = P ./ rsum





    return(P)
end
