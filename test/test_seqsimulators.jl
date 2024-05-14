using Test
using BioSequences  # Assuming SamplerWeighted and nucleotide types are defined here

@testset "simulate_sequence Function Tests" begin

    # Test 1: Basic functionality
    seq = simulate_sequence(100, [0.25, 0.25, 0.25, 0.25])
    @test length(seq) == 100
    @test all(in.(seq, Ref([DNA_A, DNA_C, DNA_G, DNA_T])))

    # Test 2: Input validation
    @test_throws ArgumentError simulate_sequence(-1, [0.25, 0.25, 0.25, 0.25])
    @test_throws ArgumentError simulate_sequence(100, [0.25, 0.25, 0.25])  # Only three frequencies
    @test_throws ArgumentError simulate_sequence(100, [0.1, 0.2, 0.3, 0.5])  # Sum not equal to 1
    @test_throws ArgumentError simulate_sequence(100, [-0.1, 0.5, 0.3, 0.3])  # Negative frequency

    # Test 3: Distribution test
    large_seq = simulate_sequence(10000, [0.1, 0.2, 0.3, 0.4])
    counts = countmap(large_seq)
    @test isapprox(counts[DNA_A] / 10000, 0.1, atol=0.05)
    @test isapprox(counts[DNA_C] / 10000, 0.2, atol=0.05)
    @test isapprox(counts[DNA_G] / 10000, 0.3, atol=0.05)
    @test isapprox(counts[DNA_T] / 10000, 0.4, atol=0.05)

    # Test 4: Edge cases
    edge_seq = simulate_sequence(1, [1, 0, 0, 0])
    @test edge_seq[1] == DNA_A
end


@testset "assign_rates Function Tests" begin

    # Test 1: Basic functionality with standard inputs
    model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0.1, substitution_model=JC())
    rates, weights = assign_rates(100, model)
    @test typeof(rates) <: Array
    @test typeof(weights) <: Dict

    # Test 2: Input validation
    @test_throws ErrorException assign_rates(100, SiteModel(mutation_rate=-0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0.1, substitution_model=JC()))
    @test_throws ErrorException assign_rates(100, SiteModel(mutation_rate=0.1, gamma_shape=-1.0, gamma_category_count=4, proportion_invariant=0.1, substitution_model=JC()))
    @test_throws ArgumentError assign_rates(100, SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=-1, proportion_invariant=0.1, substitution_model=JC()))
    @test_throws ErrorException assign_rates(100, SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=-0.1, substitution_model=JC()))

    # Test 3: Zero gamma category count
    zero_gamma_model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=0, proportion_invariant=0.0, substitution_model=JC())
    rates, weights = assign_rates(100, zero_gamma_model)
    @test all(x -> x == 0.1, rates)

    # Test 4: Gamma distribution integration
    gamma_model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0.1, substitution_model=JC())
    rates, _ = assign_rates(100, gamma_model)
    @test length(unique(rates)) > 1  # Expect more than one unique rate due to gamma distribution

end


using Test
# Assuming your DNA types and related functions are defined in a module, e.g., BioSequences
# using BioSequences 

@testset "update_site Function Tests" begin
    # Setup for tests
    weights = Dict(
        DNA_A => [0.1, 0.2, 0.3, 0.4],
        DNA_C => [0.4, 0.3, 0.2, 0.1],
        DNA_G => [0.25, 0.25, 0.25, 0.25],
        DNA_T => [0.3, 0.3, 0.2, 0.2]
    )

    # Test 1: Basic functionality
    @test update_site(DNA_A, weights) in nucleotides

    # Test 2: Correct handling of transition probabilities
    # Perform a large number of simulations to statistically test the distribution
    results = Dict(nuc => 0 for nuc in nucleotides)
    for _ in 1:10000
        result = update_site(DNA_G, weights)
        results[result] += 1
    end
    @test isapprox(results[DNA_A], 2500, atol=300)
    @test isapprox(results[DNA_C], 2500, atol=300)
    @test isapprox(results[DNA_G], 2500, atol=300)
    @test isapprox(results[DNA_T], 2500, atol=300)

    # Test 3: Handling weights that do not sum to exactly 1 (should still work if close enough)
    almost_correct_weights = copy(weights)
    almost_correct_weights[DNA_A] = [0.1, 0.2, 0.3, 0.3999]
    @test update_site(DNA_A, almost_correct_weights) in nucleotides

    # Test 6: Edge cases with zero probabilities
    zero_prob_weights = Dict(DNA_A => [1.0, 0.0, 0.0, 0.0])
    @test update_site(DNA_A, zero_prob_weights) == DNA_A
end


@testset "compute_transition_weights! Tests" begin
    # Setup for typical use case
    model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0.1, substitution_model=HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4]))
    Q = rate_matrix(model)
    λ, V, V⁻¹ = decompose(Q)
    P = zeros(4, 4)
    weights = Dict{DNA, Vector{Float64}}(DNA_A => zeros(4), DNA_C => zeros(4), DNA_G => zeros(4), DNA_T => zeros(4))
    μ = 0.1
    Δt = 1.0

    # Test basic functionality
    @testset "Basic Functionality" begin
        compute_transition_weights!(P, weights, μ, Δt, λ, V, V⁻¹)
        @test all(x -> sum(x) ≈ 1.0, values(weights))  # Ensure probabilities sum to 1
        @test all(x ≥ 0 for x in P)  # Ensure no negative probabilities
    end

    # Edge Cases
    @testset "Edge Cases" begin
        # Mutation rate or time interval is zero
        compute_transition_weights!(P, weights, 0.0, Δt, λ, V, V⁻¹)
        @test isapprox(P, I, atol=1e-5)

        compute_transition_weights!(P, weights, μ, 0.0, λ, V, V⁻¹)
        @test isapprox(P, I, atol=1e-5)
        # @test all(isapprox.([d for d in diag(P)], 1, atol=1e-5))  # Identity matrix expected when Δt is zero
    end

    # Numerical Stability
    @testset "Numerical Stability" begin
        λ_large = 1e6 * λ  # Large eigenvalues might cause instability
        compute_transition_weights!(P, weights, μ, Δt, λ_large, V, V⁻¹)
        @test all(isfinite, P)  # No Inf or NaN values
    end

    # Error Handling
    @testset "Error Handling" begin
        @test_throws DimensionMismatch compute_transition_weights!(P, weights, μ, Δt, λ[1:3], V, V⁻¹)  # Mismatched dimensions
    end
end


@testset "simulate_sequences! Function Tests" begin
    tree = rand(Nonultrametric(10))  # Assuming existence of a helper function to create trees
    site_model = SiteModel(mutation_rate=1.0, gamma_category_count=4, gamma_shape=0.5, proportion_invariant=0.1, substitution_model=JC())

    @testset "Tree Simulation" begin
        simulate_sequences!(tree, 100, site_model)
        for node in traversal(tree, preorder)
            @test length(node.data["sequence"]) == 100
            @test all(nuc -> nuc in nucleotides, node.data["sequence"])
        end
    end

    @testset "Error Handling" begin
        empty_tree = RootedTree()
        @test_throws ArgumentError simulate_sequences!(empty_tree, 100, site_model)
        @test_throws ArgumentError simulate_sequences!(tree, -1, site_model)
    end
end



@testset "propagate_sequence Function tests" begin

    # Define the different SiteModel configurations
    constant_rate_model = SiteModel(mutation_rate=0.1, gamma_shape=0., gamma_category_count=0, proportion_invariant=0., substitution_model=HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4]))
    invariant_sites_model = SiteModel(mutation_rate=0.1, gamma_shape=0., gamma_category_count=0, proportion_invariant=0.4, substitution_model=HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4]))
    gamma_variation_model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0., substitution_model=HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4]))
    combined_model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0.1, substitution_model=HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4]))

    # Helper function to compute weights and rates
    function prepare_simulation_parameters(model, n)
        rates, weights = assign_rates(n, model)
        Q = rate_matrix(model.substitution_model)
        λ, V, V⁻¹ = decompose(Q)
        P = zeros(4, 4)
        return rates, weights, λ, V, V⁻¹, P
    end

    # Test lengths for sequences
    seq_length = 100  # Length of DNA sequence for testing

    # Define a sample DNA sequence
    seq = sample([DNA_A, DNA_C, DNA_G, DNA_T], seq_length, replace=true)  # Assuming DNA_X are defined types

    @testset "Constant Mutation Rate Tests" begin

        rates, weights, λ, V, V⁻¹, P = prepare_simulation_parameters(constant_rate_model, seq_length)

        @testset "Basic Evolution" begin
            evolved_seq = propagate_sequence(seq, rates, 1.0, λ, V, V⁻¹, P, weights)
            @test length(evolved_seq) == length(seq)
            @test evolved_seq != seq  # Ensure mutation has occurred
        end
    
        @testset "Zero Time Evolution" begin
            evolved_seq = propagate_sequence(seq, rates, 0.0, λ, V, V⁻¹, P, weights)
            @test evolved_seq == seq  # No evolution when Δt is zero
        end
    
        @testset "Very Small Mutation Rate" begin
            small_rate_model = SiteModel(mutation_rate=1e-5, gamma_shape=0., gamma_category_count=0, proportion_invariant=0., substitution_model=HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4]))
            rates, weights, λ, V, V⁻¹, P = prepare_simulation_parameters(small_rate_model, seq_length)
            evolved_seq = propagate_sequence(seq, rates, 1.0, λ, V, V⁻¹, P, weights)
            @test evolved_seq == seq  # Expect no practical change due to extremely low mutation rate
        end
    
        @testset "Stress Test for Large Sequences" begin
            large_seq = repeat(seq, 10)  # Increase sequence length
            evolved_seq = propagate_sequence(large_seq, rates, 1.0, λ, V, V⁻¹, P, weights)
            @test length(evolved_seq) == length(large_seq)
        end
    end
    

    @testset "Invariant Sites Tests" begin

        seq = repeat(seq, 10)
        rates, weights, λ, V, V⁻¹, P = prepare_simulation_parameters(invariant_sites_model, length(seq))

        @testset "Proportion of Invariant Sites" begin
            evolved_seq = propagate_sequence(seq, rates, 10.0, λ, V, V⁻¹, P, weights)
            # Determine the proportion of unchanged nucleotides
            unchanged_count = sum(evolved_seq .== seq)
            @test unchanged_count / length(seq) ≥ invariant_sites_model.proportion_invariant
        end
    
        @testset "Zero Proportion Invariant" begin
            no_invariant_model = SiteModel(mutation_rate=0.1, gamma_shape=0., gamma_category_count=0, proportion_invariant=0., substitution_model=HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4]))
            rates, weights, λ, V, V⁻¹, P = prepare_simulation_parameters(no_invariant_model, seq_length)
            evolved_seq = propagate_sequence(seq, rates, 1.0, λ, V, V⁻¹, P, weights)
            @test evolved_seq != seq  # Expect changes across the sequence
        end
    
        @testset "Stress Test with Large Sequences" begin
            large_seq = repeat(seq, 10)  # Increase sequence length
            evolved_seq = propagate_sequence(large_seq, rates, 10.0, λ, V, V⁻¹, P, weights)
            @test length(evolved_seq) == length(large_seq)
            # Check if the invariant sites have not changed
            unchanged_large_count = sum(evolved_seq .== large_seq)
            @test unchanged_large_count / length(large_seq) ≥ invariant_sites_model.proportion_invariant
        end
    end
    

    @testset "Gamma Variation Tests" begin

        gamma_variation_model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0., substitution_model=HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4]))
        seq_length = 100
        sequence = sample([DNA_A, DNA_C, DNA_G, DNA_T], seq_length, replace=true)
        rates, weights, λ, V, V⁻¹, P = prepare_simulation_parameters(gamma_variation_model, seq_length)


        @testset "Rate Distribution Validation" begin
            # Verify that rates follow the expected gamma distribution
            variable_rates = filter(r -> r[2] != 0, rates)
            # Check the mean and variance of the sampled rates
            expected_mean = gamma_variation_model.mutation_rate
            expected_variance = gamma_variation_model.mutation_rate^2 / gamma_variation_model.gamma_shape
            @test isapprox(mean([r[2] for r in variable_rates]), expected_mean, atol=0.1)
            @test isapprox(var([r[2] for r in variable_rates]), expected_variance, atol=0.1)
        end
    
        @testset "Evolution with Gamma Distributed Rates" begin
            evolved_seq = propagate_sequence(sequence, rates, 10.0, λ, V, V⁻¹, P, weights)
            @test length(evolved_seq) == length(sequence)
            @test evolved_seq != seq  # Ensure mutation has occurred
        end
    end
    

    @testset "Combined Invariant Sites and Gamma Variation Tests" begin

        model_combined = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0.1, substitution_model=HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4]))
        seq_length = 100
        sequence = sample([DNA_A, DNA_C, DNA_G, DNA_T], seq_length, replace=true)
        rates, weights, λ, V, V⁻¹, P = prepare_simulation_parameters(model_combined, seq_length)


        @testset "Correct Handling of Invariant and Variable Sites" begin
            evolved_seq = propagate_sequence(sequence, rates, 1.0, λ, V, V⁻¹, P, weights)
            @test length(evolved_seq) == length(sequence)
            # Verify that the proportion of unchanged sites is at least as large as the proportion invariant
            unchanged_count = sum(evolved_seq .== sequence)
            @test unchanged_count / length(sequence) ≥ model_combined.proportion_invariant
        end
    
        @testset "Verification of Gamma Distributed Rates" begin
            # Filter out invariant sites and verify the distribution of the mutation rates on variable sites
            variable_rates = filter(r -> r[2] != 0, rates)
            expected_mean = model_combined.mutation_rate / (1.0 - model_combined.proportion_invariant)
            expected_variance = expected_mean^2 / model_combined.gamma_shape
            @test isapprox(mean([r[2] for r in variable_rates]), expected_mean, atol=0.1)
            @test isapprox(var([r[2] for r in variable_rates]), expected_variance, atol=0.1)
        end
    end
    

end