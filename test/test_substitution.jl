using LinearAlgebra

# Test for rate_matrix function
@testset "rate_matrix Function Tests" begin
    # Test with valid inputs
    π = [0.25, 0.25, 0.25, 0.25]
    R = [0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0]
    Q = SeqSim.rate_matrix(π, R)
    @test size(Q) == (4, 4)
    @test all(isapprox.(sum(Q, dims=1), 0.; atol=1e-10)) # Columns of Q should sum to 0
    @test issymmetric(Q * Diagonal(π))  # Q should be time-reversible
    @test all(isapprox.(Q * π, 0.0; atol=1e-10))  # Q should annihilate π

    # Test with non-uniform π
    π = [0.1, 0.2, 0.3, 0.4]
    Q = SeqSim.rate_matrix(π, R)
    @test size(Q) == (4, 4)
    @test all(isapprox.(sum(Q, dims=1), 0.; atol=1e-10)) # Columns of Q should sum to 0
    @test issymmetric(Q * Diagonal(π))  # Q should be time-reversible
    @test all(isapprox.(Q * π, 0.0; atol=1e-10))  # Q should annihilate π

    # Test with edge case: uniform R
    R = ones(4, 4) - I
    Q = SeqSim.rate_matrix(π, R)
    @test size(Q) == (4, 4)
    @test all(isapprox.(sum(Q, dims=1), 0.; atol=1e-10)) # Columns of Q should sum to 0
    @test issymmetric(Q * Diagonal(π))  # Q should be time-reversible
    @test all(isapprox.(Q * π, 0.0; atol=1e-10))  # Q should annihilate π

    # Test with edge case: diagonal R
    R = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    Q = SeqSim.rate_matrix(π, R)
    @test size(Q) == (4, 4)
    @test isapprox(Q, zeros(4, 4); atol=1e-10)  # Q should be zero matrix
end


# Test for JC model
@testset "JC Model Tests" begin
    # Test default constructor
    jc_model = SeqSim.JC()
    @test jc_model.π == fill(0.25, 4)
    @test size(jc_model.Q) == (4, 4)
    @test all(isapprox.(sum(jc_model.Q, dims=1), 0.; atol=1e-10))  # Columns of Q should sum to 0
    @test length(jc_model.λ) == 4
    @test size(jc_model.V) == (4, 4)
    @test size(jc_model.V⁻¹) == (4, 4)
    @test isapprox(jc_model.V * Diagonal(jc_model.λ) * jc_model.V⁻¹, jc_model.Q; atol=1e-10)  # Eigen decomposition check

    # Test numerical stability of decomposition
    poorly_conditioned_matrix = [1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0]
    λ, V, V⁻¹ = SeqSim.decompose(poorly_conditioned_matrix)
    @test length(λ) == 4
    @test size(V) == (4, 4)
    @test size(V⁻¹) == (4, 4)
end

# Test for F81 model
@testset "F81 Model Tests" begin
    # Test constructor with valid input
    f81_model = SeqSim.F81([0.1, 0.2, 0.3, 0.4])
    @test f81_model.π == [0.1, 0.2, 0.3, 0.4]
    @test size(f81_model.Q) == (4, 4)
    @test length(f81_model.λ) == 4
    @test size(f81_model.V) == (4, 4)
    @test size(f81_model.V⁻¹) == (4, 4)
    @test isapprox(f81_model.V * Diagonal(f81_model.λ) * f81_model.V⁻¹, f81_model.Q; atol=1e-10)  # Eigen decomposition check

    # Test constructor with invalid input
    @test_throws ArgumentError SeqSim.F81([-0.1, 0.2, 0.3, 0.6])  # Negative frequency
    @test_throws ArgumentError SeqSim.F81([0.1, 0.2, 0.3])        # Incorrect length
    @test_throws ArgumentError SeqSim.F81([0.1, 0.2, 0.3, 0.5])   # Frequencies don't sum to 1

    # Test numerical stability of decomposition
    poorly_conditioned_matrix = [1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0]
    λ, V, V⁻¹ = SeqSim.decompose(poorly_conditioned_matrix)
    @test length(λ) == 4
    @test size(V) == (4, 4)
    @test size(V⁻¹) == (4, 4)
end

# Test for K2P model
@testset "K2P Model Tests" begin
    # Test constructor with valid input
    k2p_model = SeqSim.K2P(2.0)
    @test k2p_model.π == fill(0.25, 4)
    @test k2p_model.κ == 2.0
    @test size(k2p_model.Q) == (4, 4)
    @test length(k2p_model.λ) == 4
    @test size(k2p_model.V) == (4, 4)
    @test size(k2p_model.V⁻¹) == (4, 4)
    @test isapprox(k2p_model.V * Diagonal(k2p_model.λ) * k2p_model.V⁻¹, k2p_model.Q; atol=1e-10)  # Eigen decomposition check

    # Test constructor with invalid input
    @test_throws ArgumentError SeqSim.K2P(-1.0)  # Negative κ
    @test_throws ArgumentError SeqSim.K2P(0.0)   # Zero κ

    # Test numerical stability of decomposition
    poorly_conditioned_matrix = [1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0]
    λ, V, V⁻¹ = SeqSim.decompose(poorly_conditioned_matrix)
    @test length(λ) == 4
    @test size(V) == (4, 4)
    @test size(V⁻¹) == (4, 4)
end


# Test for HKY model
@testset "HKY Model Tests" begin
    # Test constructor with valid input
    hky_model = SeqSim.HKY([0.1, 0.2, 0.3, 0.4], 2.0)
    @test hky_model.π == [0.1, 0.2, 0.3, 0.4]
    @test hky_model.κ == 2.0
    @test size(hky_model.Q) == (4, 4)
    @test length(hky_model.λ) == 4
    @test size(hky_model.V) == (4, 4)
    @test size(hky_model.V⁻¹) == (4, 4)
    @test isapprox(hky_model.V * Diagonal(hky_model.λ) * hky_model.V⁻¹, hky_model.Q; atol=1e-10)  # Eigen decomposition check

    # Test constructor with invalid input
    @test_throws ArgumentError SeqSim.HKY([-0.1, 0.2, 0.3, 0.6], 2.0)  # Negative frequency
    @test_throws ArgumentError SeqSim.HKY([0.1, 0.2, 0.3], 2.0)        # Incorrect length
    @test_throws ArgumentError SeqSim.HKY([0.1, 0.2, 0.3, 0.5], 2.0)   # Frequencies don't sum to 1
    @test_throws ArgumentError SeqSim.HKY([0.1, 0.2, 0.3, 0.4], -1.0)  # Negative κ
    @test_throws ArgumentError SeqSim.HKY([0.1, 0.2, 0.3, 0.4], 0.0)   # Zero κ

    # Test numerical stability of decomposition
    poorly_conditioned_matrix = [1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0]
    λ, V, V⁻¹ = SeqSim.decompose(poorly_conditioned_matrix)
    @test length(λ) == 4
    @test size(V) == (4, 4)
    @test size(V⁻¹) == (4, 4)
end


# Test for GTR model
@testset "GTR Model Tests" begin
    # Test constructor with valid input
    rates = [0 1 2 1; 1 0 1 2; 2 1 0 1; 1 2 1 0]
    π = [0.25, 0.25, 0.25, 0.25]
    gtr_model = SeqSim.GTR(π, rates)
    @test gtr_model.π == π
    @test gtr_model.rates == rates
    @test size(gtr_model.Q) == (4, 4)
    @test isapprox(sum(gtr_model.Q, dims=2), zeros(4); atol=1e-10)  # Rows of Q should sum to 0
    @test length(gtr_model.λ) == 4
    @test size(gtr_model.V) == (4, 4)
    @test size(gtr_model.V⁻¹) == (4, 4)
    @test isapprox(gtr_model.V * Diagonal(gtr_model.λ) * gtr_model.V⁻¹, gtr_model.Q; atol=1e-10)  # Eigen decomposition check

    # Test constructor with invalid input
    invalid_rates = [0 1 2; 1 0 1; 2 1 0]  # Not 4x4
    @test_throws ArgumentError SeqSim.GTR(π, invalid_rates)  # Invalid rates matrix size
    invalid_rates = [0 1 2 1; 1 0 1 2; 2 1 0 1; 1 3 1 0.5]  # Not time-reversible
    @test_throws ArgumentError SeqSim.GTR(π, invalid_rates)  # Non-time-reversible rates matrix
    invalid_π = [0.1, 0.2, 0.3]  # Incorrect length
    @test_throws ArgumentError SeqSim.GTR(invalid_π, rates)  # Invalid π length
    invalid_π = [-0.1, 0.3, 0.4, 0.4]  # Negative frequency
    @test_throws ArgumentError SeqSim.GTR(invalid_π, rates)  # Negative π values
    invalid_π = [0.1, 0.2, 0.3, 0.5]  # Frequencies don't sum to 1
    @test_throws ArgumentError SeqSim.GTR(invalid_π, rates)  # Invalid π sum

    # Test numerical stability of decomposition
    poorly_conditioned_matrix = [1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0]
    λ, V, V⁻¹ = SeqSim.decompose(poorly_conditioned_matrix)
    @test length(λ) == 4
    @test size(V) == (4, 4)
    @test size(V⁻¹) == (4, 4)
end
