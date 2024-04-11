export logistic

function logistic(
    supremum::Vector{Float64}, 
    steepness::Float64
    )

    midpoint = supremum .* 0.5


    ## g is the logistic transformation ([-inf, inf] => [0.0, supremum])
    g(x::Vector{T}) where {T <: Real} = begin
        supremum ./ (1 .+ exp.((-1) .* steepness .* (x .- midpoint)))        
    end

    ## h is the inverse operation of g ([0.0, supremum] => [-inf, inf])
    h(y::Vector{T}) where {T <: Real} = begin
        midpoint .+ (log.(y) .- log.(supremum .- y)) ./ steepness
    end
    
    return(g, h) ## returns two functions
end

function logistic(
        lower::Vector{Float64},
        upper::Vector{Float64},
        steepness::Float64
    )

    midpoint = (lower .+ upper) .* 0.5
    ## cheating a bit here, so the inverse function doesn't go to +/- Inf
    #ϵ = 1e-20 ## an offset

    g(x::Vector{T}) where {T <: Real} = begin
        (upper .- lower) ./ (1 .+ exp.(-steepness .* (x .- midpoint))) .+ lower
    end

    h(y::Vector{T}) where {T <: Real} = begin
        #y = [maximum([yi, li + ϵ]) for (yi, li) in zip(y, lower)]
        #y = [minimum([yi, li - ϵ]) for (yi, li) in zip(y, upper)]
        midpoint .- log.((upper .- lower) ./ (y .- lower) .- 1) ./ steepness
    end

    return(g,h)
end

