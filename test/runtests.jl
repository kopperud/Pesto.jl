using Pesto
using Test

@testset "Likelihood regression" begin

    λ = [0.2, 0.1]
    μ = [0.05, 0.1]
    η = 0.005
    model = BDSconstant(λ, μ, η)

    include("../docs/src/primates.jl");


    logl = logL_root(model, primates)

    # Write your tests here.
    threshold = 0.001
    @test abs(logl - -693.1624874593532) < threshold

end

@testset "Branch rates" begin
    λ = [0.2, 0.1]
    μ = [0.05, 0.1]
    η = 0.005
    model = BDSconstant(λ, μ, η)

    include("../docs/src/primates.jl");

    rates = tree_rates(primates, model);

    target = [0.1088643102855757, 0.12872942956813035, 0.18662616176687286, 0.19813663825163044, 0.19705776696087068, 0.19766996131326536, 0.19710110266036726, 0.19710110266036726, 0.19864070515647547, 0.1984029272470074, 0.19712588682309887, 0.19712588682309887, 0.19892770622412512, 0.19896287964727813, 0.19764486196519887, 0.19879042162272886, 0.19854426120380406, 0.19809736853583723, 0.19809736853583723, 0.19808507061531294, 0.1987173751049797, 0.1974617912799596, 0.1974617912799596, 0.17868850526251295, 0.17998699690565034, 0.1802814904233706, 0.18097822664971505, 0.18069381717127436, 0.18069381717127436, 0.1795041005146961]


    threshold = 0.0001
    @test sum((rates[1:30,:mean_lambda] .- target).^2) < threshold

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
