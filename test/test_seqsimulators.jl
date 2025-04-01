using Test
using StaticArrays

# Test suite for `update_transition_weights!`
@testset "update_transition_weights!" begin
    using LinearAlgebra

    substitution_model = SeqSim.HKY([0.3, 0.2, 0.3, 0.2], 2.0)
    Q = substitution_model.Q
    λ = substitution_model.λ
    V = substitution_model.V
    V⁻¹ = substitution_model.V⁻¹

    # Test 1: Basic functionality with a single category
    transition_weights = [zeros(4, 4)]
    Δt = 0.1
    μ = [0.5]

    expected = V * Diagonal(exp.(μ[1] * λ * Δt)) * V⁻¹
    SeqSim.update_transition_weights!(transition_weights, Δt, μ, λ, V, V⁻¹)
    @test transition_weights[1] ≈ expected

    # Test 2: Multiple categories
    transition_weights = [zeros(4, 4) for _ in 1:3]
    Δt = 0.2
    μ = [0.3, 0.6, 0.9]

    for cat in 1:3
        expected = V * Diagonal(exp.(μ[cat] * λ * Δt)) * V⁻¹
        SeqSim.update_transition_weights!(transition_weights, Δt, μ, λ, V, V⁻¹)
        @test transition_weights[cat] ≈ expected
    end

    # Test 3: Edge case with zero mutation rates
    transition_weights = [zeros(4, 4)]
    Δt = 0.5
    μ = [0.0]

    expected = V * Diagonal(exp.(μ[1] * λ * Δt)) * V⁻¹
    SeqSim.update_transition_weights!(transition_weights, Δt, μ, λ, V, V⁻¹)
    @test transition_weights[1] ≈ expected

    # Test 4: Edge case with zero time interval
    transition_weights = [zeros(4, 4)]
    Δt = 0.0
    μ = [0.5]

    expected = V * Diagonal(exp.(μ[1] * λ * Δt)) * V⁻¹
    SeqSim.update_transition_weights!(transition_weights, Δt, μ, λ, V, V⁻¹)
    @test transition_weights[1] ≈ expected
end


# Test suite for `update_site`
@testset "update_site" begin
    # Test 1: Basic functionality with uniform probabilities
    weights = [0.25 0.25 0.25 0.25;
               0.25 0.25 0.25 0.25;
               0.25 0.25 0.25 0.25;
               0.25 0.25 0.25 0.25]
    for nucl_in in UInt8(1):UInt8(4)
        results = [SeqSim.update_site(nucl_in, weights) for _ in 1:10^4]
        counts = countmap(results)
        for nucl_out in UInt8(1):UInt8(4)
            @test abs(counts[nucl_out] / 10^4 - 0.25) ≤ 0.05
        end
    end

    # Test 2: Transition probabilities favoring one nucleotide
    weights = [0.7 0.1 0.1 0.1;
               0.1 0.7 0.1 0.1;
               0.1 0.1 0.7 0.1;
               0.1 0.1 0.1 0.7]
    for nucl_in in UInt8(1):UInt8(4)
        results = [SeqSim.update_site(nucl_in, weights) for _ in 1:10^4]
        counts = countmap(results)
        @test abs(counts[nucl_in] / 10^4 - 0.7) ≤ 0.05
        for nucl_out in UInt8(1):UInt8(4)
            if nucl_out != nucl_in
                @test abs(counts[nucl_out] / 10^4 - 0.1) ≤ 0.05
            end
        end
    end

    # Test 3: Edge case with deterministic transitions
    weights = [1.0 0.0 0.0 0.0;
               0.0 1.0 0.0 0.0;
               0.0 0.0 1.0 0.0;
               0.0 0.0 0.0 1.0]
    for nucl_in in UInt8(1):UInt8(4)
        results = [SeqSim.update_site(nucl_in, weights) for _ in 1:10^4]
        @test all(result == nucl_in for result in results)
    end

    # Test 4: Edge case with zero probabilities (invalid input)
    weights = [0.0 0.0 0.0 0.0;
               0.0 0.0 0.0 0.0;
               0.0 0.0 0.0 0.0;
               0.0 0.0 0.0 0.0]
    for nucl_in in UInt8(1):UInt8(4)
        @test_throws BoundsError SeqSim.update_site(nucl_in, weights)
    end

    # Test 5: Non-uniform probabilities
    weights = [0.4 0.3 0.2 0.1;
               0.3 0.4 0.1 0.2;
               0.2 0.1 0.4 0.3;
               0.1 0.2 0.3 0.4]
    for nucl_in in UInt8(1):UInt8(4)
        results = [SeqSim.update_site(nucl_in, weights) for _ in 1:10^4]
        counts = countmap(results)
        for nucl_out in UInt8(1):UInt8(4)
            expected_prob = weights[nucl_out, nucl_in]
            @test abs(counts[nucl_out] / 10^4 - expected_prob) ≤ 0.05
        end
    end
end


@testset "update_sequence!" begin
    # Test 1: Basic functionality with a single category
    seq_in = UInt8[1, 2, 3, 4]
    Δt = 0.1
    μ = [0.5]
    variable_sites = [[1, 2, 3, 4]]
    substitution_model = SeqSim.HKY([0.3, 0.2, 0.3, 0.2], 2.0)
    Q = substitution_model.Q
    λ = substitution_model.λ
    V = substitution_model.V
    V⁻¹ = substitution_model.V⁻¹
    transition_weights = [zeros(4, 4)]

    expected = SeqSim.update_transition_weights!(transition_weights, Δt, μ, λ, V, V⁻¹)
    seq_out = SeqSim.update_sequence!(transition_weights, seq_in, Δt, μ, variable_sites, λ, V, V⁻¹)
    @test length(seq_out) == length(seq_in)
    @test all(nucl in UInt8(1):UInt8(4) for nucl in seq_out)

    # Test 2: Multiple categories
    seq_in = UInt8[1, 2, 3, 4, 1, 2, 3, 4]
    Δt = 0.2
    μ = [0.3, 0.6]
    variable_sites = [[1, 2, 3, 4], [5, 6, 7, 8]]
    transition_weights = [zeros(4, 4) for _ in 1:2]

    expected = SeqSim.update_transition_weights!(transition_weights, Δt, μ, λ, V, V⁻¹)
    seq_out = SeqSim.update_sequence!(transition_weights, seq_in, Δt, μ, variable_sites, λ, V, V⁻¹)
    @test length(seq_out) == length(seq_in)
    @test all(nucl in UInt8(1):UInt8(4) for nucl in seq_out)

    # Test 3: Edge case with zero mutation rates
    seq_in = UInt8[1, 2, 3, 4]
    Δt = 0.5
    μ = [0.0]
    variable_sites = [[1, 2, 3, 4]]
    transition_weights = [zeros(4, 4)]

    seq_out = SeqSim.update_sequence!(transition_weights, seq_in, Δt, μ, variable_sites, λ, V, V⁻¹)
    @test seq_out == seq_in

    # Test 4: Edge case with zero time interval
    seq_in = UInt8[1, 2, 3, 4]
    Δt = 0.0
    μ = [0.5]
    variable_sites = [[1, 2, 3, 4]]
    transition_weights = [zeros(4, 4)]

    seq_out = SeqSim.update_sequence!(transition_weights, seq_in, Δt, μ, variable_sites, λ, V, V⁻¹)
    @test seq_out == seq_in

    # Test 5: Non-uniform variable sites
    seq_in = UInt8[1, 2, 3, 4, 1, 2, 3, 4]
    Δt = 0.3
    μ = [0.4, 0.7]
    variable_sites = [[1, 3, 5, 7], [2, 4, 6, 8]]
    transition_weights = [zeros(4, 4) for _ in 1:2]

    seq_out = SeqSim.update_sequence!(transition_weights, seq_in, Δt, μ, variable_sites, λ, V, V⁻¹)
    @test length(seq_out) == length(seq_in)
    @test all(nucl in UInt8(1):UInt8(4) for nucl in seq_out)

    # Test 6: Edge case with empty variable sites
    seq_in = UInt8[1, 2, 3, 4]
    Δt = 0.1
    μ = Vector{Float64}()
    variable_sites = Vector{Vector{Int64}}()
    transition_weights = Vector{Matrix{Float64}}()

    seq_out = SeqSim.update_sequence!(transition_weights, seq_in, Δt, μ, variable_sites, λ, V, V⁻¹)
    @test seq_out == seq_in
end
