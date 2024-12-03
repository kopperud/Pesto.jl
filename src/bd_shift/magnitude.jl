export magnitude

@doc raw"""
    magnitude(model, data)

Computes the overall point estimate of the magnitude of the rate shifts in the tree, for a given model. The equation for the magnitude is
```math
\text{mag} = \frac{1}{\sum_j \sum_i \hat{N}_{ij}} \sum_j \sum_i (r_i - r_j) \hat{N}_{ij},
```
where ``\hat{N}_{ij}`` is the number of estimated diversification rate shifts (summed across all branches) that depart from the rate category `j` and arrive in the rate category `i` (in the direction of old to young).

Example:

```julia
phy = readtree(Pesto.path("primates.tre"))
sampling_probability = 0.635
primates = SSEdata(phy, sampling_probability)

λ = [0.05, 0.15, 0.25]
μ = [0.03, 0.08, 0.05]
η = 0.003

model = BDSconstant(λ, μ, η)

mag = magnitude(model, primates)
println(mag)
```
"""
function magnitude(model::Model, data::SSEdata)
    N = state_shifts(model, data)
    mag = magnitude(model, N)

    return(mag)
end

function magnitude(model::Model, N::Array{Float64, 3})
    N = sum(N, dims = 1)[1,:,:]

    mag = magnitude(model, N)
    return(mag)
end

function magnitude(model::Model, N::Array{Float64, 2})
    r = model.λ .- model.μ
    Δr = r .- r' 

    nom = sum(Δr .* N)
    denom = sum(N)
    mag = nom / denom

    return(mag)
end

