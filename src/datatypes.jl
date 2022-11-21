export SSEconstant, SSEdata, BDconstant

struct SSEconstant <: Distributions.ContinuousUnivariateDistribution
    λ
    μ
    η
end

struct SSEdata
    state_space
    trait_data
    edges
    tiplab
    node_depth
    ρ
    branch_lengths
    branching_times
    po
end

struct BDconstant <: Distributions.ContinuousUnivariateDistribution
    λ
    μ
end

struct SSEresult
    phy
    lambda
    mu
end
