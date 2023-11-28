export plot_preorder!

function plot_preorder!(
    h::Array{Float64,2}, 
    #v, 
    t::Float64, 
    i::Vector{Int64}, 
    r::Vector{Float64}, 
    data::SSEdata, 
    edge_index::Int64
    )

    node_index = data.edges[edge_index,2]
    bl = data.branch_lengths[edge_index]
    h[edge_index,1:2] = [t, t + bl]
    #n_species = length(data.tiplab)

    # if is internal
    if node_index > data.Nnode+1
        left, right = findall(data.edges[:,1] .== node_index)
        plot_preorder!(h, t+bl, i, r, data, left)
        plot_preorder!(h, t+bl, i, r, data, right)
        h[edge_index,3] = (h[left,3] + h[right,3]) / 2
    # if is tip 
    else
        h[edge_index,3] = i[1]
        i[1] += 1
    end
    h[edge_index,4] = r[edge_index]
end

export coordinates

function coordinates(data::SSEdata, r::Vector{Float64})
    n_species = length(data.tiplab)
    root_index = n_species+1
    left, right = findall(data.edges[:,1] .== root_index)
    n_edges = length(data.branch_lengths)
    h = zeros(n_edges, 4)

#    v = Vector{Float64}[]
    i = [1]

    t = 0.0
    plot_preorder!(h, t, i, r, data, left)
    plot_preorder!(h, t, i, r, data, right)

    #h = h[sortperm(h[:, 1]), :]
    return(h)
end

export plot_tree

#RecipesBase.@userplot plot

RecipesBase.@recipe function f(data::SSEdata)
    r = zeros(size(data.edges)[1])
    h = coordinates(data, r)
    #h = h[sortperm(h[:, 1]), :]
    h = sortslices(h, dims = 1, by = x -> (x[1],x[2]))
    

    for i in 1:size(h)[1]
        x = h[i,1:2]
        y = [h[i,3], h[i,3]]
        RecipesBase.@series begin
            seriestype := :line
            color := :black
            label := ""
            grid := false
            yaxis := false
            xaxis := false
            x, y
        end
        #plot!(p, , , label = "", )

        ## vertical branches        
        if i % 2 > 0
            x = [h[i,1], h[i+1,1]]
            y = [h[i,3], h[i+1,3]]

            RecipesBase.@series begin
                seriestype := :line
                color := :black
                label := ""
                x, y
            #plot!(p, x, y, label = "", color = :black)
            end
        end

        ## tip labels

    end
end

import ColorSchemes 

RecipesBase.@recipe function f(
    data::SSEdata, 
    df::DataFrames.DataFrame,
    palette = ColorSchemes.linear_ternary_blue_0_44_c57_n256;
    #palette = ColorSchemes.linear_ternary_red_0_50_c52_n256;
    rate = "mean_lambda"
    )

    df2 = DataFrames.sort(df, :edge)
    r = df2[!,Symbol(rate)][2:end]


    #framestyle := [:none, :none]

    #legend := false
    #layout := RecipesBase.@layout [a{0.95w} b]
    #layout := (2,1)
    
    
    h = coordinates(data, r)
    h = sortslices(h, dims = 1, by = x -> (x[1],x[2]))

    colors = ColorSchemes.get(palette, h[:,4], :extrema)
    
    for i in 1:size(h)[1]
        x = h[i,1:2]
        y = [h[i,3], h[i,3]]
        RecipesBase.@series begin
            seriestype := :line
            color := :black
            label := ""
            grid := false
            yaxis := false
            xaxis := false
            color := colors[i]
            title := rate
            subplot := 1
            x, y
        end
        #plot!(p, , , label = "", )

        ## vertical branches    
        if i % 2 > 0 ## only for internal branches

            bottom = h[i+1,3]
            top = h[i,3]

            mid = (bottom+top) / 2

            ## branch going up
            x = [h[i,1], h[i,1]]
            y = [mid, top]

            RecipesBase.@series begin
                seriestype := :line
                color := :black
                label := ""
                color := colors[i]
                subplot := 1
                x, y
            end

            ## branch going down
            y = [mid, bottom]

            RecipesBase.@series begin
                seriestype := :line
                color := :black
                label := ""
                color := colors[i+1]
                subplot := 1
                x, y
            end
        end
        ## tip labels
    end
#= 
    cmap = :thermal
    RecipesBase.@series begin
        seriestype := :heatmap
        clims := extrema(r)
        #subplot := 2
        framestyle := :none
        c := cmap
        cbar := true
        #lims := (-1,0)
    end =#
end