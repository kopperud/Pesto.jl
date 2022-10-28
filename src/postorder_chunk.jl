export postorder_chunk

function postorder_chunk(data)
    visited_nodes = Set()
    po = []

    global edges = data.edges
    ntips = length(data.tiplab)
    root_node = length(data.tiplab)+1

    append!(po, [collect(1:ntips)])

    while !(root_node in visited_nodes)
        nrows = size(edges)[1]
        is_descendant = [edges[j,2] in po[end] for j in 1:nrows]

        if length(is_descendant) == 0
            break
        end

        push!(visited_nodes, edges[is_descendant,2]...)

        ## check nodes have both their descendants already traversed
        po1 = Int64[]
        for j in edges[:,2]
            children = data.edges[data.edges[:,1] .== j,2]
            if all([child in visited_nodes for child in children])
                if !(j in visited_nodes)
                    append!(po1, j)
                end
            end
        end
        if length(po1) > 0 
            append!(po, [po1])
        end

        global edges = edges[.!is_descendant, :]
    end 

    return(po)
end


function levelorder(data)
    edges = data.edges
    node = length(data.tiplab)+1
    visits = Int64[]

    q = DataStructures.Queue{Int}()
    DataStructures.enqueue!(q, node)

    while !DataStructures.isempty(q)
        node = DataStructures.dequeue!(q)

        append!(visits, node)
        children = edges[edges[:,1] .== node,2]

        l = length(children)
        if l > 0
            for child in children
                DataStructures.enqueue!(q, child)
            end
        end
    end
    return(visits)
end

