## read newick

function readfile(filename)
    io = open(filename, "r")
    s = read(io, String)
    close(io)

    return(s)
end

function readnewick(filename)
    s = readfile(filename)

    y = [item[1] for item in findall(";", s)]

    s = replace(s, r"^_+|_+$" => "")
    s = replace(s, r"[ \t]" => "")

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

function findsplit(tokens)
    global ps = 0
    for (i, token) in enumerate(tokens)
        if token == "("
            ps += 1
        elseif token == ")"
            ps -= 1
        end
        
        if (token == ",") & (ps == 0)
            return(i)
        end
    end
    throw("split not found") 
end

function partition(tokens)
    comma = findsplit(tokens)

    left = tokens[1:comma-1]
    right = tokens[1+comma:end]

    return (left, right)
end

function internaledge!(edges, tokens, idx, node)
    tokens = tokens[2:end-2]

    edges[idx[1],1] = node
    edges[idx[1],2] = maximum(edges)+1
    node = edges[idx[1],2]
    idx[1] += 1    
    left, right = partition(tokens)

    for branch in [left, right]
        if !isempty(branch)
            if branch[end][1] != ':'
                println(idx)
                terminaledge!(edges, branch[1], idx, node)
            else
                internaledge!(edges, branch, idx, node)
            end
        end
    end
end

function terminaledge!(edges, s, idx, node)
    edges[idx[1],1] = node
    edges[idx[1],2] = idx[2]

    idx[1] += 1
    idx[2] += 1
end