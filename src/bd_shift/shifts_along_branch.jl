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
        model::BhDhModel, 
        data::SSEdata;
        summarize = :geometric,
        condition = [:mrca, :survival],
        height = maximum(data.node_depth),
        tspan = range(height, 0.0; length = 500),
    )

    Ds, Fs = backwards_forwards_pass(model, data; condition = condition);

    y = Float64[]
    times = Float64[]

    nbranches = size(data.edges)[1]
    oldest_times = [Ds[edge_idx].t[end] for edge_idx in 1:nbranches]
    youngest_times = [Ds[edge_idx].t[1] for edge_idx in 1:nbranches]

    K = length(model.λ)
    r = -(model.η/(K-1.0))

    #for t in range(height,0.0; length = 500)
    for t in tspan
        active_branch_indices = which_active_branches(t, youngest_times, oldest_times)
        n_active = length(active_branch_indices)

        ΔN_dts = Float64[]
        for edge_idx in active_branch_indices
            Dt = Ds[edge_idx](t)[:,2]
            Ft = Fs[edge_idx](t)
            St = ancestral_state_probability(Dt, Ft, t)
            dN_dt = - r * (sum(Dt .* sum(St ./ Dt)) -1)

            if dN_dt > 0.0
                push!(ΔN_dts, dN_dt)
            else
                push!(ΔN_dts, 1e-8) 
            end
        end
        
        if !isempty(ΔN_dts)
            if summarize == :arithmetic
                # arithmetic mean
                mean_ΔN_dt = sum(ΔN_dts) / length(ΔN_dts)
            elseif summarize == :geometric
                # geometric mean
                mean_ΔN_dt = exp(sum(log.(ΔN_dts) / length(ΔN_dts)))
            end

            append!(y, mean_ΔN_dt)
            append!(times, t) 
        end
    end

    return(times, y)
end
