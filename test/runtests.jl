using Pesto
using Test

@testset "Pesto.jl" begin

    λ = [0.2, 0.1]
    μ = [0.05, 0.1]
    η = 0.005
    model = SSEconstant(λ, μ, η)

    include("../docs/src/primates.jl");


    logl = logL_root(model, primates)

    # Write your tests here.
    @test abs(logl - -693.1624874593532) < 0.001

end
