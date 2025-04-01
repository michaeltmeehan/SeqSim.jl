using Test
using Random
using Distributions

@testset "assign_rates" begin
    # Test 1: Basic functionality with gamma categories
    @testset "Basic functionality with gamma categories" begin
        sequence_length = 1000
        proportion_invariant = 0.2
        mutation_rate = 1e-3
        gamma_shape = 0.5
        gamma_category_count = 4

        variable_sites, μ = SeqSim.assign_rates(sequence_length, proportion_invariant, mutation_rate, gamma_shape, gamma_category_count)

        # Check the number of invariant and variable sites
        num_invariant_sites = Int(floor(proportion_invariant * sequence_length))
        num_variable_sites = sequence_length - num_invariant_sites
        @test length(reduce(vcat, variable_sites)) == num_variable_sites
        @test length(variable_sites) == gamma_category_count

        # Check mutation rates
        @test length(μ) == gamma_category_count
        @test all(rate > 0 for rate in μ)
        @test length.(variable_sites) ⋅ μ ./ sequence_length ≈ mutation_rate
    end

    # Test 2: No invariant sites
    @testset "No invariant sites" begin
        sequence_length = 1000
        proportion_invariant = 0.0
        mutation_rate = 1e-3
        gamma_shape = 1.0
        gamma_category_count = 3

        variable_sites, μ = SeqSim.assign_rates(sequence_length, proportion_invariant, mutation_rate, gamma_shape, gamma_category_count)

        @test length(reduce(vcat, variable_sites)) == sequence_length
        @test length(variable_sites) == gamma_category_count
        @test length(μ) == gamma_category_count
        @test isapprox(length.(variable_sites) ⋅ μ ./ sequence_length, mutation_rate, rtol=1e-3)
    end

    # Test 3: All invariant sites
    @testset "All invariant sites" begin
        sequence_length = 1000
        proportion_invariant = 1.0
        mutation_rate = 1e-3
        gamma_shape = 1.0
        gamma_category_count = 3

        variable_sites, μ = SeqSim.assign_rates(sequence_length, proportion_invariant, mutation_rate, gamma_shape, gamma_category_count)

        @test length(reduce(vcat, variable_sites)) == 0
        @test length(variable_sites) == gamma_category_count
        @test isempty(reduce(vcat, variable_sites))
        @test isnan(μ[1])
    end

    # Test 4: Single gamma category
    @testset "Single gamma category" begin
        sequence_length = 1000
        proportion_invariant = 0.3
        mutation_rate = 1e-3
        gamma_shape = 2.0
        gamma_category_count = 1

        variable_sites, μ = SeqSim.assign_rates(sequence_length, proportion_invariant, mutation_rate, gamma_shape, gamma_category_count)

        num_variable_sites = sequence_length - Int(floor(proportion_invariant * sequence_length))
        @test length(variable_sites[1]) == num_variable_sites
        @test length(μ) == 1
        @test length.(variable_sites) ⋅ μ ./ sequence_length ≈ mutation_rate
    end

    # Test 5: Edge case with zero gamma categories
    @testset "Zero gamma categories" begin
        sequence_length = 1000
        proportion_invariant = 0.3
        mutation_rate = 1e-3
        gamma_shape = 2.0
        gamma_category_count = 0

        variable_sites, μ = SeqSim.assign_rates(sequence_length, proportion_invariant, mutation_rate, gamma_shape, gamma_category_count)

        num_variable_sites = sequence_length - Int(floor(proportion_invariant * sequence_length))
        @test length(variable_sites[1]) == num_variable_sites
        @test length(μ) == 1
    end

    # Test 6: Check that all sites are accounted for
    @testset "All sites accounted for" begin
        sequence_length = 1000
        proportion_invariant = 0.25
        mutation_rate = 1e-3
        gamma_shape = 1.0
        gamma_category_count = 5

        variable_sites, μ = SeqSim.assign_rates(sequence_length, proportion_invariant, mutation_rate, gamma_shape, gamma_category_count)

        all_variable_sites = reduce(vcat, variable_sites)
        @test unique(all_variable_sites) == all_variable_sites
        @test length(all_variable_sites) == sequence_length - Int(floor(proportion_invariant * sequence_length))
    end
end


@testset "SiteModel Tests" begin
    # Mock substitution model for testing
    struct MockSubstitutionModel <: SubstitutionModel end
    mock_model = MockSubstitutionModel()

    # Test 1: Basic functionality
    @testset "Basic functionality" begin
        sequence_length = 1000
        mutation_rate = 1e-3
        gamma_category_count = 4
        gamma_shape = 0.5
        proportion_invariant = 0.2

        site_model = SiteModel(sequence_length, mutation_rate, gamma_category_count, gamma_shape, proportion_invariant, mock_model)

        @test site_model.sequence_length == sequence_length
        @test site_model.mutation_rate == mutation_rate
        @test site_model.gamma_category_count == gamma_category_count
        @test site_model.gamma_shape == gamma_shape
        @test site_model.proportion_invariant == proportion_invariant
        @test site_model.substitution_model == mock_model
        @test length(site_model.variable_sites) == gamma_category_count
        @test length(site_model.μ) == gamma_category_count
        @test all(rate > 0 for rate in site_model.μ)
    end

    # Test 2: Invalid inputs
    @testset "Invalid inputs" begin
        @test_throws ArgumentError SiteModel(-100, 1e-3, 4, 0.5, 0.2, mock_model)
        @test_throws ArgumentError SiteModel(100, -1e-3, 4, 0.5, 0.2, mock_model)
        @test_throws ArgumentError SiteModel(100, 1e-3, -4, 0.5, 0.2, mock_model)
        @test_throws ArgumentError SiteModel(100, 1e-3, 4, -0.5, 0.2, mock_model)
        @test_throws ArgumentError SiteModel(100, 1e-3, 4, 0.5, -0.2, mock_model)
        @test_throws ArgumentError SiteModel(100, 1e-3, 4, 0.5, 1.2, mock_model)
    end
end
