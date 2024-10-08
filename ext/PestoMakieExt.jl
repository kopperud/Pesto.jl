module PestoMakieExt

import Makie, Pesto, DataFrames
import ColorSchemes

## MAKIE recipe

function Pesto.treeplot!(
    ax::Makie.Axis,
    data::Pesto.SSEdata;
    tip_label = false,
    tip_label_size = 5
)

    n_edges = size(data.edges)[1]
    n_tips = length(data.tiplab)

    r = zeros(n_edges)

    h = Pesto.coordinates(data, r)
    horizontal_lines = [
        Makie.Point(h[i,1], h[i,3]) => Makie.Point(h[i,2], h[i,3]) for i in 1:n_edges
    ]
    h = sortslices(h, dims = 1, by = x -> (x[1],x[2]))

    upward_lines = similar(horizontal_lines, 0)
    downward_lines = similar(horizontal_lines, 0)
    r_up = Float64[]
    r_down = Float64[]

    for i in 1:n_edges
        if i % 2 > 0 ## only for internal edges
            top = h[i+1,3]
            bottom = h[i,3]
            #bottom, top = extrema(h[i:i+1, 3])

            mid = (bottom+top) / 2

            ## branch going up
            start = Makie.Point(h[i,1], mid)
            finish = Makie.Point(h[i,1], top)
            append!(upward_lines, [start => finish])

            ## branch going down
            finish = Makie.Point(h[i,1], bottom)
            append!(downward_lines, [start => finish])

            append!(r_up, h[i+1,4])
            append!(r_down, h[i,4])
        end
    end

    ## actually do the plotting
    Makie.linesegments!(ax, downward_lines, color = :black, linewidth = 1)
    Makie.linesegments!(ax, upward_lines, color = :black, linewidth = 1)
    Makie.linesegments!(ax, horizontal_lines, color = :black, linewidth = 1)
        
    Makie.hidespines!(ax)
    Makie.hidedecorations!(ax)


    if tip_label
        for i in 1:n_tips
            tp = data.tiplab[i]
            Makie.text!(ax, tp, fontsize = tip_label_size, 
                    position = (maximum(h[:,2])+0.5, i), 
                    align = (:left, :center))
        end
        tip_label_x_offset = 5.0
    else
        tip_label_x_offset = 0.0
    end
    top_offset = 0.05*n_tips
    Makie.ylims!(ax, (0.0, n_tips + top_offset))
end

function Pesto.treeplot!(
    ax::Makie.Axis,
    data::Pesto.SSEdata, 
    rates::DataFrames.DataFrame,
    rate_name = "mean_netdiv",
    cmap = Makie.cgrad(:Spectral, 5, categorical = true);
    tip_label = false,
    tip_label_size = 5
)

    n_edges = size(data.edges)[1]
    n_tips = length(data.tiplab)

    @assert rate_name in names(rates)
    df2 = DataFrames.sort(rates, :edge)
    r = df2[!,Symbol(rate_name)][2:end]

    h = Pesto.coordinates(data, r)
    horizontal_lines = [
        Makie.Point(h[i,1], h[i,3]) => Makie.Point(h[i,2], h[i,3]) for i in 1:n_edges
    ]
    h = sortslices(h, dims = 1, by = x -> (x[1],x[2]))

    upward_lines = similar(horizontal_lines, 0)
    downward_lines = similar(horizontal_lines, 0)
    r_up = Float64[]
    r_down = Float64[]

    for i in 1:n_edges
        if i % 2 > 0 ## only for internal edges
            top = h[i+1,3]
            bottom = h[i,3]
            #bottom, top = extrema(h[i:i+1, 3])

            mid = (bottom+top) / 2

            ## branch going up
            start = Makie.Point(h[i,1], mid)
            finish = Makie.Point(h[i,1], top)
            append!(upward_lines, [start => finish])

            ## branch going down
            finish = Makie.Point(h[i,1], bottom)
            append!(downward_lines, [start => finish])

            append!(r_up, h[i+1,4])
            append!(r_down, h[i,4])
        end
    end

    ## actually do the plotting
    r_extrema = extrema(r)
    ## do the color bar
    #if !isnothing(rate_name)
    #cmap = Makie.cgrad(:Spectral, 5, categorical = true)
    

    Makie.linesegments!(ax, downward_lines, color = r_down, colormap = cmap, colorrange = r_extrema, linewidth = 1)
    Makie.linesegments!(ax, upward_lines, color = r_up, colormap = cmap, colorrange = r_extrema, linewidth = 1)
    Makie.linesegments!(ax, horizontal_lines, color = r, colormap = cmap, linewidth = 1)

    
    Makie.hidespines!(ax)
    Makie.hidedecorations!(ax)

    if tip_label
        for i in 1:n_tips
            tp = data.tiplab[i]
            Makie.text!(ax, tp, fontsize = tip_label_size, 
                    position = (maximum(h[:,2])+0.5, i), 
                    align = (:left, :center))
        end
        tip_label_x_offset = 5.0
    else
        tip_label_x_offset = 0.0
    end
    ## 5% buffer
    top_offset = 0.05*n_tips
    Makie.ylims!(ax, (0.0, n_tips + top_offset))
end

function Pesto.treeplot(
    data::Pesto.SSEdata, 
    rates::DataFrames.DataFrame,
    rate_name::String = "mean_netdiv";
    cmap = Makie.cgrad(:Spectral, 5, categorical = true),
    tip_label = false,
    tip_label_size = 5
)
    if !isnothing(rate_name)
        @assert rate_name in names(rates) "rate name must be equal to one of the column names in the rates data frame"
    end

    rates = deepcopy(rates)
    sort!(rates, :edge)

    #cmap = Makie.cgrad(:Spectral, 5, categorical = true)

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1])
    Pesto.treeplot!(ax, data, rates, rate_name, cmap; tip_label = tip_label, tip_label_size = tip_label_size)
    if !isnothing(rate_name)
        Makie.Colorbar(fig[0,1], limits = extrema(rates[2:end, Symbol(rate_name)]), label = rate_name, colormap = cmap, vertical = false)
        Makie.rowgap!(fig.layout, 2.0)
    end
    return(fig)
end

function Pesto.treeplot(
    data::Pesto.SSEdata;
    tip_label = false,
    tip_label_size = 5
)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1])
    Pesto.treeplot!(ax, data; tip_label = tip_label, tip_label_size = tip_label_size)
    return(fig)
end


##########################################
#
#          plot logl surface
#          
##########################################


function Pesto.plot_loglsurface(
        vectors::Vector{Vector{Float64}},
        logls::Vector{Matrix{Float64}};
        offsets = [10.0, 10.0, 10.0]
    )
    logls = deepcopy(logls)

    λs, μs, ηs = vectors

    for (i, offset) in enumerate(offsets)
        max = maximum(logls[i])    
        #logls[i][:,:] .= logls[i][:,:] .- max ## set top to equal 0

        logls[i][(logls[i] .< (max-offset))] .= max -offset
    end

    fig = Makie.Figure(size = (300, 600))

    viewmode = :stretch
    elevation = 0.3 * pi
    #            (left,right, bottom?, top?)
    protrusions = (50, 20, 0, 0)

    ax1 = Makie.Axis3(fig[1,1], xlabel = "μ", ylabel = "η", zlabel = "logl"; viewmode, elevation, protrusions)
    ax2 = Makie.Axis3(fig[2,1], xlabel = "λ", ylabel = "η", zlabel = "logl"; viewmode, elevation, protrusions) 
    ax3 = Makie.Axis3(fig[3,1], xlabel = "λ", ylabel = "μ", zlabel = "logl"; viewmode, elevation, protrusions) 

    methods = [
               Makie.surface!, 
               Makie.contour3d!
              ]
    for method! in methods 
        method!(ax1, μs, ηs, logls[1])
        method!(ax2, λs, ηs, logls[2])
        method!(ax3, λs, μs, logls[3])
    end

    ## draw plane with a mesh (two triangles)
    # (x, y, z)
    ma = maximum(logls[3]) 
    mi = minimum(logls[3]) 

    p1 = Makie.Point3(minimum(λs), minimum(λs), ma)
    p2 = Makie.Point3(maximum(λs), maximum(λs), ma)
    p3 = Makie.Point3(minimum(λs), minimum(λs), mi)
    p4 = Makie.Point3(maximum(λs), maximum(λs), mi)
    
    Makie.mesh!(ax3, [p1, p2, p3], color = (:red, 0.3), shading = Makie.automatic, transprent = true)
    Makie.mesh!(ax3, [p2, p3, p4], color = (:red, 0.3), shading = Makie.automatic, transprent = true)

    Makie.rowgap!(fig.layout, 0.0)
    Makie.colgap!(fig.layout, 0.0)

    return(fig)
end


end
