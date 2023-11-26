export plot_preorder!

function plot_preorder!(
    h::Array{Float64,2}, 
    #v, 
    t::Float64, 
    i::Vector{Int64}, 
    data::SSEdata, 
    edge_index::Int64
    )

    node_index = data.edges[edge_index,2]
    bl = data.branch_lengths[edge_index]
    h[edge_index,1:2] = [t, t + bl]
    n_species = length(data.tiplab)

    # if is internal
    if node_index > data.Nnode+1
        left, right = findall(data.edges[:,1] .== node_index)
        plot_preorder!(h, t+bl, i, data, left)
        plot_preorder!(h, t+bl, i, data, right)
        h[edge_index,3] = (h[left,3] + h[right,3]) / 2
    # if is tip 
    else
        h[edge_index,3] = i[1]
        i[1] += 1
    end
end

export coordinates

function coordinates(data::SSEdata)
    n_species = length(data.tiplab)
    root_index = n_species+1
    left, right = findall(data.edges[:,1] .== root_index)
    n_edges = length(data.branch_lengths)
    h = zeros(n_edges, 3)

#    v = Vector{Float64}[]
    i = [1]

    t = 0.0
    plot_preorder!(h, t, i, data, left)
    plot_preorder!(h, t, i, data, right)

    #h = h[sortperm(h[:, 1]), :]
    return(h)
end

export plot_tree

#RecipesBase.@userplot plot

RecipesBase.@recipe function f(data::SSEdata)
    h = coordinates(data)
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