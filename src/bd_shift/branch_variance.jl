export branch_variance
export branch_mean

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

export brvar

function brvar(
        model,
        data,
        branch_index,
        num_time_steps,
    )
    Ds, Fs = backwards_forwards_pass(model, data)
    Ss = ancestral_state_probabilities(data, Ds, Fs)
    E = extinction_probability(model, data)

    D = Ds[branch_index]
    F = Fs[branch_index]
    S = Ss[branch_index]

    tmax = maximum(D.t)
    tmin = minimum(D.t)

    times = collect(range(tmax, tmin; length = num_time_steps))
    Δt = times[2] - times[1]

    Dt = D(tmax)[:,2]
    Ft = F(tmax)[:,2]
    St = ancestral_state_probability(Dt, Ft, tmax)

    means = Float64[]
    vars = Float64[]
    prob = Float64[]
    time_index = 1
    EX = sum(S(times[time_index]) .* model.λ)

    for state_index in 1:length(model.λ)
        m = model.λ[state_index] / num_time_steps
        v = (EX - model.λ[state_index])^2 / num_time_steps
        p = S(tmax)[state_index]
        brvar!(model, D, E, S, Δt, times, time_index, state_index, num_time_steps, means, vars, prob, m, v, p)
    end

    return(means, vars, prob)
end

function brvar!(
        model::BhDhModel, 
        D, 
        E, 
        S,
        Δt::Float64, 
        times::Vector{Float64}, 
        time_index::Int64, 
        state_index::Int64,
        num_time_steps::Int64, 
        means::Vector{Float64}, 
        vars::Vector{Float64}, 
        prob::Vector{Float64}, 
        m::Float64,
        v::Float64,
        p::Float64)
    t = times[time_index]
    P = Pmatrix(model, D, E, t, Δt)
    EX = sum(S(times[time_index+1]) .* model.λ)

    for i in 1:length(model.λ)
        m1 = m + model.λ[i] / num_time_steps
        v1 = v + (EX - model.λ[i])^2 / num_time_steps
        p1 = p * P[i,state_index] 

        if time_index < (num_time_steps - 1)
            brvar!(model, D, E, S, Δt, times, time_index+1, i, num_time_steps, means, vars, prob, m1, v1, p1)
        else
            push!(means, m1)
            push!(vars, v1)
            push!(prob, p1)
        end
    end
end

function branch_variance(model, data, branch_index, n)
    Ds, Fs = backwards_forwards_pass(model, data)
    Ss = ancestral_state_probabilities(data, Ds, Fs)
    E = extinction_probability(model, data)

    D = Ds[branch_index]
    F = Fs[branch_index]
    E = extinction_probability(model, data)

    domain = Ds[branch_index].t
    tmin = minimum(domain)
    tmax = maximum(domain)
    n_categories = length(model.λ)

    times = collect(range(tmax, tmin; length = n))
    
    var = 0.0
    Δt = times[2] - times[1]

    for time_idx in 1:(n-1)
        t = times[time_idx]
        Dt = D(t)[:,2]
        Ft = F(t)[:,2]
        St = ancestral_state_probability(Dt, Ft, t)

        ## Transition probability matrix, given state j
        P = Pmatrix(model, D, E, t, Δt)

        ## 1st moment
        ## E[X(t+Δt)]
        m1 = sum(sum(model.λ[j] * St[j] * P[i,j] for i in 1:4) for j in 1:4)
        ## 2nd moment
        ## E[X(t+Δt)^2]
        m2 = sum(sum(model.λ[j] * model.λ[j] * St[j] * P[i,j] for i in 1:4) for j in 1:4)
        ## variance of time step t 
        ## (assuming independence, i.e. 0 covariance)
        v = m2 - m1^2

        ## add to the total
        var += v / (n^2)
    end

    return(var)
end

export foobaz

function foobaz(model, data, branch_index)
    Ds, Fs = backwards_forwards_pass(model, data)
    Ss = ancestral_state_probabilities(data, Ds, Fs)
    E = extinction_probability(model, data)

    D = Ds[branch_index]
    F = Fs[branch_index]
    
    tmax = maximum(D.t)
    tmin = minimum(D.t)

    t = tmax
    Δt = -5.0
    Dt = D(t)[:,2]
    Ft = F(t)[:,2]

    times = collect(range(tmax, tmin; length = 100))
    Δt = times[2] - times[1]
    println(string("Δt: ", Δt))

    cv1s = Float64[]
    cv2s = Float64[]
    for i in 1:99
        t = times[i]
        P = Pmatrix(model, D, E, t, Δt)
        St = ancestral_state_probability(Dt, Ft, t)

        cv1 = sum(sum(model.λ[i] * model.λ[j] * St[j] * P[i,j] for i in 1:4) for j in 1:4)
        cv2 = sum(sum(model.λ[i] * model.λ[j] * Ss[branch_index](t)[j] * Ss[branch_index](t+Δt)[i] for i in 1:4) for j in 1:4)

        push!(cv1s, cv1)
        push!(cv2s, cv2)
    end
        
    return(cv1s, cv2s)
end

export Amatrix, Pmatrix

function Amatrix(model, E, K, t)
    Q = -LinearAlgebra.I(K) .* model.η .+ (1 .- LinearAlgebra.I(K)) .* (model.η/(K-1))
    A = LinearAlgebra.diagm(- (model.λ .+ model.μ) .+ 2 .* model.λ .* E(t)) .+ Q
    return(A)
end

function Pmatrix(model, D, E, t, Δt)
    K = length(model.λ)
    A = Amatrix(model, E, K, t)
    Dt = D(t)[:,2]

    P_unnorm = (LinearAlgebra.I(K) .- Δt .* A) .* (Dt * ones(K)')
    csum = ones(K) * sum(P_unnorm, dims = 1) ## column sum
    P = P_unnorm ./ csum
    return(P)
end
