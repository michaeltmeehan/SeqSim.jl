function empirical_state_frequencies(states::Vector{UInt8})
    return [count(==(UInt8(i)), states) / length(states) for i in 1:4]
end


@testset "JC propagation matches transition probabilities" begin
    site_model = SiteModel(MersenneTwister(101), 20_000, 0.7, 0, 0.0, 0.0, JC())
    propagator = SequencePropagator(site_model)
    root = fill(UInt8(1), site_model.sequence_length)
    Δt = 0.8

    child = propagator(MersenneTwister(102), root, Δt)
    SeqSim.update_transition_weights!(
        propagator.transition_weights,
        Δt,
        site_model.μ,
        propagator.decomposition.λ,
        propagator.decomposition.V,
        propagator.decomposition.V⁻¹,
    )

    expected = vec(propagator.transition_weights[1][:, 1])
    observed = empirical_state_frequencies(child)
    @test all(abs.(observed .- expected) .< 0.02)
end


@testset "long-time F81 propagation approaches stationary frequencies" begin
    π = [0.1, 0.2, 0.3, 0.4]
    site_model = SiteModel(MersenneTwister(201), 30_000, 1.0, 0, 0.0, 0.0, F81(π))
    propagator = SequencePropagator(site_model)
    root = fill(UInt8(1), site_model.sequence_length)

    child = propagator(MersenneTwister(202), root, 12.0)
    observed = empirical_state_frequencies(child)
    @test all(abs.(observed .- π) .< 0.02)
end


@testset "small elapsed time produces rare changes under low rate" begin
    site_model = SiteModel(MersenneTwister(301), 10_000, 0.01, 0, 0.0, 0.0, JC())
    propagator = SequencePropagator(site_model)
    root = repeat(UInt8[1, 2, 3, 4], div(site_model.sequence_length, 4))

    child = propagator(MersenneTwister(302), root, 0.01)
    changed = count(root .!= child)
    @test changed <= 5
end


@testset "invariant-site mixture preserves exactly the intended sites" begin
    site_model = SiteModel(MersenneTwister(401), 200, 4.0, 3, 0.8, 0.25, HKY([0.25, 0.25, 0.25, 0.25], 2.0))
    propagator = SequencePropagator(site_model)
    root = repeat(UInt8[1, 2, 3, 4], div(site_model.sequence_length, 4))

    child = propagator(MersenneTwister(402), root, 2.0)
    variable_sites = sort(reduce(vcat, site_model.variable_sites))
    invariant_sites = setdiff(1:site_model.sequence_length, variable_sites)

    @test length(invariant_sites) == floor(Int, site_model.proportion_invariant * site_model.sequence_length)
    @test child[invariant_sites] == root[invariant_sites]
    @test count(root[variable_sites] .!= child[variable_sites]) > 0
end


@testset "explicit RNG controls site assignment reproducibly" begin
    args = (101, 0.2, 4, 0.75, 0.2, JC())
    model_a = SiteModel(MersenneTwister(501), args...)
    model_b = SiteModel(MersenneTwister(501), args...)
    model_c = SiteModel(MersenneTwister(502), args...)

    @test model_a.variable_sites == model_b.variable_sites
    @test model_a.μ == model_b.μ
    @test model_a.variable_sites != model_c.variable_sites

    sizes = length.(model_a.variable_sites)
    @test maximum(sizes) - minimum(sizes) <= 1
    @test issorted(model_a.μ)
    @test all(diff(model_a.μ) .> 0.0)
    @test isapprox(sizes ⋅ model_a.μ / model_a.sequence_length, model_a.mutation_rate; rtol=0.02)
end
