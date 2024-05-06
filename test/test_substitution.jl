@testset "Substitution Models Tests" begin
    # Test for Jukes-Cantor Model
    @testset "Jukes-Cantor Model" begin
        model = JC()
        @test model.π == fill(0.25, 4)  # Test default base frequencies
        rm = rate_matrix(model)
        @test size(rm) == (4, 4)  # Ensure the matrix is of correct size
        @test all(isapprox.(sum(rm, dims=1), 0.; atol=1e-5))  # Sum of each row should be zero
        @test all(isapprox.(rm * model.π, 0; atol=1e-5))    # Stationarity
        @test issymmetric(rm * diagm(model.π))   # Detailed balance
    end

    # Test for Felsenstein 81 Model
    @testset "Felsenstein 81 Model" begin
        model = F81(π=[0.1, 0.2, 0.3, 0.4])
        @test model.π == [0.1, 0.2, 0.3, 0.4]  # Check if π values are set correctly
        rm = rate_matrix(model)
        @test size(rm) == (4, 4)  # Ensure the matrix is of correct size
        @test all(isapprox.(sum(rm, dims=1), 0.; atol=1e-5))  # Check for row sum to zero
        @test rm[2, 1] == 0.2  # Rates should match target frequencies
        @test all(isapprox.(rm * model.π, 0; atol=1e-5))    # Stationarity
        @test issymmetric(rm * diagm(model.π))   # Detailed balance
    end

    # Test for Kimura 2-Parameter Model
    @testset "Kimura 2-Parameter Model" begin
        κ = 2.0
        model = K2P(κ=κ, π = fill(0.25, 4))
        @test model.κ == κ  # Check if κ is set correctly
        rm = rate_matrix(model)
        @test size(rm) == (4, 4)  # Ensure the matrix is of correct size
        @test all(isapprox.(sum(rm, dims=1), 0.; atol=1e-5)) # Row sum should be zero
        @test rm[3, 1] == κ * rm[2, 1]  # Transitions should be κ times more likely than transversions
        @test all(isapprox.(rm * model.π, 0; atol=1e-5))    # Stationarity
        @test issymmetric(rm * diagm(model.π))  # Detailed balance
    end

    # Test for Hasegawa-Kishino-Yano Model
    @testset "Hasegawa-Kishino-Yano Model" begin
        κ = 2.0
        model = HKY(κ=κ, π=[0.1, 0.2, 0.3, 0.4])
        @test model.κ == κ && model.π == [0.1, 0.2, 0.3, 0.4]  # Check κ and π are set correctly
        rm = rate_matrix(model)
        @test size(rm) == (4, 4)  # Ensure the matrix is of correct size
        @test all(isapprox.(sum(rm, dims=1), 0.; atol=1e-5))  # Ensure the row sums to zero
        @test rm[3, 1] == κ * 0.3  # Transition rates check
        @test all(isapprox.(rm * model.π, 0; atol=1e-5))    # Stationarity
        @test issymmetric(rm * diagm(model.π))   # Detailed balance
    end

    # Test for General Time Reversible Model
    # @testset "General Time Reversible Model" begin
    #     rate_matrix_input = [0 1 2 1; 1 0 1 2; 2 1 0 1; 1 2 1 0]
    #     model = GTR(rate_matrix=rate_matrix_input, π=[0.25, 0.25, 0.25, 0.25])
    #     @test model.π == [0.25, 0.25, 0.25, 0.25]  # Check base frequencies
    #     rm = rate_matrix(model)
    #     @test size(rm) == (4, 4)  # Ensure the matrix is of correct size
    #     @test all(isapprox.(sum(rm, dims=1), 0.; atol=1e-5)) # Ensure that each row sums to zero
    #     # @test isapprox(rm[1, 2], π[2])  # Test if rate matrix follows time-reversibility
        # @test all(isapprox.(rm * model.π, 0; atol=1e-5))    # Stationarity
        # @test issymmetric(rm * diagm(model.π))   # Detailed balance
    # end
end
