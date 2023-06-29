## read newick

function readfile(filename)
    io = open(filename, "r")
    s = read(io, String)
    close(io)

    return(s)
end

function readnewick(filename)
    s = readfile(filename)

    #y = [item[1] for item in findall(";", s)]
    s = replace(s, r"^_+|_+$" => "")
    s = replace(s, r"[ \t]" => "")
    s = s[findfirst('(', s):end]
    s = s[1:findfirst(';', s)]
    s = stripcomments(s)
    tokens = tokenize(s)

    ntip = sum(1 for token in tokens if token == "(") +1
    nroot = ntip+1
    idx = [1,1]
    edges = zeros(Int64, 2*(ntip-1), 2)
    edges[1,1] = nroot
    el = zeros(size(edges)[1])
    tiplabs = String[]

    tokens2 = tokens[2:end-2]
    left, right = partition(tokens2)


    for side in [left, right]
        if length(side) == 1
            terminaledge!(edges, el, tiplabs, side, idx, nroot)
        else
            internaledge!(edges, el, tiplabs, side, idx, nroot)
        end
    end

    return(edges, el, tiplabs)
end

function stripcomments(s)
    res = replace(s, r"\[.*?\]" => "")
    return(res)
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

function parse_tiplab(s)
    res = split(s, ':')[1]
    return(res)
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

function internaledge!(edges, el, tiplabs, tokens, idx, node)
    l = parse_brlen(tokens[end])
    el[idx[1]] = l
    tokens = tokens[2:end-2]

    edges[idx[1],1] = node
    edges[idx[1],2] = maximum(edges)+1
    node = edges[idx[1],2]
    idx[1] += 1
    left, right = partition(tokens)

    for branch in [left, right]
        if !isempty(branch)
            if branch[end][1] != ':'
                #println(idx)
                terminaledge!(edges, el, tiplabs, branch[1], idx, node)
            else
                internaledge!(edges, el, tiplabs, branch, idx, node)
            end
        end
    end
end

function terminaledge!(edges, el, tiplabs, s, idx, node)
    edges[idx[1],1] = node
    edges[idx[1],2] = idx[2]

    l = parse_brlen(s)
    el[idx[1]] = l

    tiplab = String(split(s, ':')[1])
    push!(tiplabs, tiplab)

    idx[1] += 1
    idx[2] += 1
end