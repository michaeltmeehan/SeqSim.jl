
@testset "mod_wrap Function Tests" begin
    # Test standard cases
    @testset "Standard Cases" begin
        @test mod_wrap(10, 3) == 1
        @test mod_wrap(18, 5) == 3
        @test mod_wrap(25, 7) == 4
        @test mod_wrap(0, 5) == 5  # Boundary case: x is zero
    end

    # Test cases where x % y == 0
    @testset "Modulo Zero Cases" begin
        @test mod_wrap(20, 5) == 5
        @test mod_wrap(12, 3) == 3
        @test mod_wrap(0, 1) == 1  # Both x and y leading to zero modulo
    end

    # Test large numbers
    @testset "Large Numbers" begin
        @test mod_wrap(100000000, 123456) == 100000000 % 123456
    end

    # Test edge cases for small and large values of y
    @testset "Edge Cases" begin
        @test mod_wrap(5, 100) == 5
        @test mod_wrap(12345, 1) == 1  # Modulo by 1
    end

    # Test error handling
    @testset "Error Handling" begin
        @test_throws ArgumentError mod_wrap(10, 0)  # Zero divisor
        @test_throws ArgumentError mod_wrap(5, -5)  # Negative divisor
    end
end


@testset "decompose Function Tests" begin
    # Test correct decomposition with a simple diagonal matrix
    @testset "Diagonal Matrix" begin
        M = diagm(0 => [1.0, 2.0, 3.0])
        λ, V, V⁻¹ = decompose(M)
        @test λ == [1.0, 2.0, 3.0]
        @test V == I  # Eigenvectors for a diagonal matrix should be the identity matrix
        @test V⁻¹ == I  # Inverse of the identity matrix is itself
    end

    # Test correct decomposition with a symmetric matrix
    @testset "Symmetric Matrix" begin
        M = [2 1; 1 2]
        λ, V, V⁻¹ = decompose(M)
        @test all(λ .≈ [1.0, 3.0])  # Eigenvalues should be approximately 1 and 3
        @test isapprox(V * diagm(0 => λ) * V⁻¹, M)  # AV = VD should hold, hence VDV⁻¹ ≈ A
    end

    # Test handling of poorly conditioned matrices
    # @testset "Poorly Conditioned Matrix" begin
    #     M = [1 1000; 0 1]  # This matrix is poorly conditioned
    #     λ, V, V⁻¹ = decompose(M)
    #     @test size(λ) == (2,)
    #     @test cond(V) > 1e12  # Condition number should be high
    #     @test V * diagm(0 => λ) * V⁻¹ ≈ M  # The decomposition should still reconstruct M approximately
    # end

    # Test error handling for non-square matrices
    @testset "Non-Square Matrix" begin
        M = [1 2 3; 4 5 6]
        @test_throws ArgumentError decompose(M)
    end
end
