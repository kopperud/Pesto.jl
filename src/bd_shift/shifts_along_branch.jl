function which_active_branches(
        t::Float64, 
        youngest_times::Vector{Float64},
        oldest_times::Vector{Float64}, 
    )

    active_lineages = (oldest_times .> t) .& (youngest_times .< t)

    active_branch_indices =  [idx for (idx, bool) in enumerate(active_lineages) if bool]
    return(active_branch_indices)
end

export shift_rate_through_time

function shift_rate_through_time(
        model::BDSconstant, 
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
