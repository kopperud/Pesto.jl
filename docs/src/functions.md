# Functions

List of available functions

```@docs
birth_death_shift(model::Model, data::SSEdata)
```

```@docs
birth_death_shift(model::Model, tree::Root)
```

```@docs
Econstant(t, λ, µ, ρ)
```

```@docs
psi(t, λ, µ, ρ)
```

```@docs
lp(λ::Vector{Float64}, μ::Vector{Float64}, data::SSEdata)
```

```@docs
estimate_constant_bdp(data::SSEdata)
```

```@docs
optimize_eta(λ::Vector{Float64}, µ::Vector{Float64}, data::SSEdata)
```

```@docs
make_descendants(data::SSEdata)
```

```@docs
make_ancestors(data::SSEdata)
```

```@docs
lrange(from::Float64, to::Float64, length::Int64)
```

```@docs
allpairwise(λ, µ)
```

```@docs
writenewick(fpath::String, data::SSEdata, rates::DataFrames.DataFrame)
```

```@docs
SSEdata(phy::Pesto.phylo, ρ::Float64)
```

```@docs
make_descendants_nodes(data::SSEdata)
```

```@docs
readtree(path::String)
```

```@docs
postorder_async(model::SSEconstant, data::SSEdata, E)
```

```@docs
postorder_sync(model::SSEconstant, data::SSEdata, E)
```

```@docs
magnitude(model::SSEconstant, data::SSEdata)
```

```@docs
pesto_twostep(tree::SSEdata)
```

```@docs
pesto_fossil(tree::Root)
```

```@docs
pesto(tree::SSEdata)
```
