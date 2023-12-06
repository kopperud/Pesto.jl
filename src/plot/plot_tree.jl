export plot_preorder!
export treeplot

function plot_preorder!(
    h::Array{Float64,2}, 
    t::Float64, 
    i::Vector{Int64}, 
    r::Vector{Float64}, 
    data::SSEdata, 
    edge_index::Int64
    )

    node_index = data.edges[edge_index,2]
    bl = data.branch_lengths[edge_index]
    h[edge_index,1:2] = [t, t + bl]

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

function coordinates(
    data::SSEdata
)
    r = zeros(size(data.edges)[1])
    coordinates(data, r)
end

function coordinates(
    data::SSEdata, 
    r::Vector{Float64}
    )

    n_species = length(data.tiplab)
    root_index = n_species+1
    left, right = findall(data.edges[:,1] .== root_index)
    n_edges = length(data.branch_lengths)
    h = zeros(n_edges, 4)

    i = [1]
    t = 0.0

    plot_preorder!(h, t, i, r, data, left)
    plot_preorder!(h, t, i, r, data, right)

    return(h)
end


## create new function
## the method is empty, but we need it initialized here,
## so that the plot extension module can create new methods in Pesto.treeplot()
function treeplot()
end

function treeplot!()
end