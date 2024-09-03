using Pesto
using Test

@testset "Likelihood regression" begin

    λ = [0.2, 0.1]
    μ = [0.05, 0.1]
    η = 0.005
    model = SSEconstant(λ, μ, η)

    include("../docs/src/primates.jl");


    logl = logL_root(model, primates)

    # Write your tests here.
    threshold = 0.001
    @test abs(logl - -693.1624874593532) < threshold

end

@testset "Tree checks" begin
    @testset "Polytomies" begin
        phy = readtree(Pesto.path("polytomy.tre"))

        Test.@test_throws PolytomyError SSEdata(phy, 1.0)
    end

    @testset "Negative branch lengths" begin
        phy = readtree(Pesto.path("negative_branch.tre"))

        Test.@test_throws NegativeBranchError SSEdata(phy, 1.0)
    end

    @testset "Non-ultrametric tree" begin
        phy = readtree(Pesto.path("nonultrametric.tre"))

        Test.@test_warn "Your tree appears to not be ultrametric. Check if this is a rounding error issue, or if it really is not ultrametric. Pesto only works for ultrametric trees." SSEdata(phy, 1.0)
    end
end
