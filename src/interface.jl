export interface

function interface(;fossil_extinction = 0.0, shifts = :equal, kwargs...)
    d = Dict{Symbol,Symbol}(kwargs)

    keys1 = [:speciation, :extinction, :rel_ext, :net_div]

    vals1 = [:zero, :constant, :varying]

    parameters = Symbol[]

    nkeys1 = 0
    for (key, val) in d
        if key in keys1
            if nkeys1 > 1
                error("only two of $keys1 can be chosen")
            end

            push!(parameters, val)
            nkeys1 += 1
        end
       
        if key == :fossilization
            nothing
        end
    end

    for v in values(d)
        if !(v in vals1)
            error("arguments must be one of $vals1")
        end
    end

    args2 = [:equal, :independent]

    for (x, y) in d
        println(x, "\t", y)
    end

    println("shifts", "\t", shifts)
end
