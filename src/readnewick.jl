## read newick

function readtree(filename)
    io = open(filename, "r")

    x = read(io, String)
    close(io)

    y = [item[1] for item in findall(";", x)]

    x = replace(x, r"^_+|_+$" => "")
    x = replace(x, r"[ \t]" => "")



    return(y)
end

function stripcomments(s)
    res = replace(s, r"\[.*?\]" => "")
end


function tokenize(s)
    tokens = String[]

    ## strip everything between square bracket
    s = stripcomments(s)
    len = length(s)

    single_tokens = Set([')', '(', ',', ';'])
    i = 1
    while i <= len
        if s[i] ∈ single_tokens
            token = string(s[i])
            append!(tokens, [token])
            i += 1
        else
            l = Int64[]
            firstcomma = findfirst(',', s[i:end])
            firstclose = findfirst(')', s[i:end])
            if !isnothing(firstcomma)
                append!(l, firstcomma-1)
            end
            if !isnothing(firstclose)
                append!(l, firstclose-1)
            end

            if !isempty(l)
                close_idx = minimum(l)
                token = s[i:close_idx+i-1]
                append!(tokens, [token])
                i += length(token)
            else
                i += 1
            end 
        end
    end
    return(tokens)
end

function parse_brlen(s)
    res = parse(Float64, split(s, ':')[end])
    return(res)
end

function buildtree(tokens)
    tr = []

    #ntip = length(findall(")", s))+1
    ntip = 1
    for token in tokens
        if token == "("
            ntip += 1
        end
    end
    node_idx = ntip + 1
    
    nbranches = ntip * 2 - 2
    edges = zeros(Int64, nbranches, 2)
    el = zeros(nbranches)
    edge_idx = 1
    tip_idx = 1
    i = 1

    addbranch(tokens, edges, el, i, edge_idx, node_idx, tip_idx)
    return(edges, el)
end

function addbranch(tokens, edges, el, i, edge_idx, node_idx, tip_idx)

    #i += 1
    single_tokens = Set([')', '(', ',', ';'])
    if tokens[i] == "("
        edges[edge_idx,:] = [node_idx, node_idx+1]
        node_idx += 1
        edge_idx += 1
        i += 1
        addbranch(tokens, edges, el, i, edge_idx, node_idx, tip_idx)
    elseif tokens[i] == ")"
        el[edge_idx] = parse_brlen(tokens[i+1])
        i += 1
#    elseif tokens[i] == ";"      
    else
        ## tip
        edges[edge_idx, :] = [node_idx, tip_idx]
        tip_idx += 1
    end
    #return(edges, el)
end

function foo()        
    if tokens[i] ∉ single_tokens
        if tokens[i][1] == ':' ## internal
            tip!(tokens, edges, el, i)
        else ## external
            addbranch!(tokens, edges, el, i)
        end
    end

    internal!(tokens, edges, el, i)
    el[i] = parse_brlen(tokens[i])    
end