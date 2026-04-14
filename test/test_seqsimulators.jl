@testset "root sequence simulation reproducibility and frequencies" begin
    frequencies = [0.1, 0.2, 0.3, 0.4]
    seq1 = rand_seq(MersenneTwister(42), 2000; frequencies=frequencies, taxon="root")
    seq2 = rand_seq(MersenneTwister(42), 2000; frequencies=frequencies, taxon="root")

    @test seq1.value == seq2.value
    @test length(seq1.value) == 2000
    @test all(base in "ACGT" for base in seq1.value)

    empirical = [count(==(base), seq1.value) / length(seq1.value) for base in "ACGT"]
    @test all(abs.(empirical .- frequencies) .< 0.04)

    @test_throws ArgumentError rand_seq(MersenneTwister(1), 0)
    @test_throws ArgumentError rand_seq(MersenneTwister(1), 10; frequencies=[0.5, 0.5])
    @test_throws ArgumentError rand_seq(MersenneTwister(1), 10; frequencies=[0.4, 0.3, 0.2, 0.2])
end

@testset "sequence propagation invariants" begin
    Random.seed!(21)
    site_model = SiteModel(16, 0.5, 2, 1.0, 0.25, HKY([0.25, 0.25, 0.25, 0.25], 2.0))
    propagator = SequencePropagator(site_model)
    root = SeqSim.encode("ACGTACGTACGTACGT")

    same = propagator(MersenneTwister(7), root, 0.0)
    @test same == root

    child1 = propagator(MersenneTwister(8), root, 1.0)
    child2 = propagator(MersenneTwister(8), root, 1.0)
    @test child1 == child2
    @test length(child1) == length(root)
    @test all(state in UInt8(1):UInt8(4) for state in child1)

    invariant_sites = setdiff(1:site_model.sequence_length, reduce(vcat, site_model.variable_sites))
    @test child1[invariant_sites] == root[invariant_sites]

    @test_throws ArgumentError propagator(MersenneTwister(1), UInt8[1, 2, 3], 0.1)
    @test_throws ArgumentError propagator(MersenneTwister(1), UInt8[1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3], 0.1)
    @test_throws ArgumentError propagator(MersenneTwister(1), root, -0.1)
end

@testset "transition weights are stochastic matrices" begin
    site_model = SiteModel(8, 0.2, 1, 1.0, 0.0, JC())
    propagator = SequencePropagator(site_model)
    SeqSim.update_transition_weights!(
        propagator.transition_weights,
        0.5,
        site_model.μ,
        propagator.decomposition.λ,
        propagator.decomposition.V,
        propagator.decomposition.V⁻¹,
    )

    weights = propagator.transition_weights[1]
    @test all(weights .>= -1e-12)
    @test all(isapprox.(sum(weights, dims=1), 1.0; atol=1e-10))
end
