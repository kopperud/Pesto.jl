export shifts_tiny_intervals

function shifts_tiny_intervals(
        model::SSE, 
        data::SSEdata,
        Ds, 
        Fs,
    )
    

    nshifts = nshifts_simple_whole_branch(model, data, Ds, Fs);
    nbranches = size(data.edges)[1]

    ## chop up N along a branch into tiny pieces ΔN_i = N(t_{i+1}) - N(t_i)
    ## ΔN_foo holds for each branch an interpolated function that
    ## represents ΔN_i for some tiny time interval Δt 
    ΔN_foo = Dict()
    for edge_idx in 1:nbranches
        a = Ds[edge_idx].t[end]
        b = Ds[edge_idx].t[1]
        tspan = (a,b)

        Δt_target = 0.01

        number_episodes = Int64(round((a-b)/Δt_target)+1)
        Δt = (a-b) / number_episodes

        times = collect(range(a, b; length = number_episodes+1))

        ΔNs = Float64[]
        for i in 2:(number_episodes+1)
            ΔN = nshifts[edge_idx](times[i]) - nshifts[edge_idx](times[i-1])
            append!(ΔNs, ΔN)
        end

        mid_times = (times[1:end-1] .+ times[2:end]) ./ 2

        if length(mid_times) > 1
            ## Flat() extrapolation means that if t is too old, then the oldest value of f(t) is used
            ## likewise, if t is too young, then the youngest value of f(t) is used
            f = Interpolations.linear_interpolation(reverse(mid_times), reverse(ΔNs), extrapolation_bc = Interpolations.Flat())
        else
            f = t -> ΔNs[1]
        end

        ΔN_foo[edge_idx] = f
    end

    return(ΔN_foo)
end


function which_active_branches(
        t::Float64, 
        youngest_times::Vector{Float64},
        oldest_times::Vector{Float64}, 
    )

    active_lineages = (oldest_times .> t) .& (youngest_times .< t)

    active_branch_indices =  [idx for (idx, bool) in enumerate(active_lineages) if bool]
    return(active_branch_indices)
end


function mean_ΔN_through_time(
        model::SSE, 
        data::SSEdata;
        Δt = 0.05,
        height = maximum(data.node_depth),
        tspan = range(height, 0.0; length = 500),
    )

    Ds, Fs = backwards_forwards_pass(model, data);

    y = Float64[]
    times = Float64[]

    nbranches = size(data.edges)[1]
    oldest_times = [Ds[edge_idx].t[end] for edge_idx in 1:nbranches]
    youngest_times = [Ds[edge_idx].t[1] for edge_idx in 1:nbranches]

    ΔN_foo = shifts_tiny_intervals(model, data, Ds, Fs)

    for t in range(height,0.0; length = 500)
        active_branch_indices = which_active_branches(t, youngest_times, oldest_times)
        n_active = length(active_branch_indices)

        ΔNs = Float64[]
        for edge_idx in active_branch_indices
           x = ΔN_foo[edge_idx](t) 
           append!(ΔNs, x)
        end
        
        if !isempty(ΔNs)
            #mean_ΔN = Statistics.mean(ΔNs)
            mean_ΔN = sum(ΔNs) / length(ΔNs)
            append!(y, mean_ΔN)    
            append!(times, t) 
        end
    end

    return(times, y)
end


function mean_ΔNdt_through_time(
        model::SSE, 
        data::SSEdata;
        height = maximum(data.node_depth),
        tspan = range(height, 0.0; length = 500),
    )

    Ds, Fs = backwards_forwards_pass(model, data);

    y = Float64[]
    times = Float64[]

    nbranches = size(data.edges)[1]
    oldest_times = [Ds[edge_idx].t[end] for edge_idx in 1:nbranches]
    youngest_times = [Ds[edge_idx].t[1] for edge_idx in 1:nbranches]

    ΔN_foo = shifts_tiny_intervals(model, data, Ds, Fs)

    K = length(model.λ)
    r = -(model.η/(K-1.0))

    for t in range(height,0.0; length = 500)
        active_branch_indices = which_active_branches(t, youngest_times, oldest_times)
        n_active = length(active_branch_indices)


        ΔN_dts = Float64[]
        for edge_idx in active_branch_indices
            Dt = Ds[edge_idx](t)
            Ft = Fs[edge_idx](t)
            St = ancestral_state_probability(Dt, Ft, t)
            dN_dt = - r * (sum(Dt .* sum(St ./ Dt)) -1)

            append!(ΔN_dts, dN_dt)
        end
        
        if !isempty(ΔN_dts)
            mean_ΔN_dt = sum(ΔN_dts) / length(ΔN_dts)
            append!(y, mean_ΔN_dt)    
            append!(times, t) 
        end
    end

    return(times, y)
end
