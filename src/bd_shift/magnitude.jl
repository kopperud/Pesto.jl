export magnitude

@doc raw"""
    magnitude(model, data)

Computes the overall point estimate of the magnitude of the rate shifts in the tree, for a given model.

Example:

```julia
phy = readtree(Pesto.path("primates.tre"))
ρ = 0.635
primates = SSEdata(phy, ρ)

λ = [0.05, 0.15]
μ = [0.03, 0.08]
η = 0.01

model = SSEconstant(λ, μ, η)

mag = magnitude(model, data)
println(mag)
```
"""
function magnitude(model::SSE, data::SSEdata)
    N = state_shifts(model, data)
    mag = magnitude(model, data, N)

    return(mag)
end

function magnitude(model::SSE, data::SSEdata, N::Array{Float64, 3})
    N = sum(N, dims = 1)[1,:,:]

    r = model.λ .- model.μ
    Δr = r .- r' 

    nom = sum(Δr .* N)
    denom = sum(N)
    mag = nom / denom

    return(mag)
end
