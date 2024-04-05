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