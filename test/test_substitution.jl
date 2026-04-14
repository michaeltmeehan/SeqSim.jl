@testset "substitution model invariants" begin
    models = [
        JC(),
        F81([0.1, 0.2, 0.3, 0.4]),
        K2P(2.0),
        HKY([0.2, 0.3, 0.3, 0.2], 3.0),
        GTR([0.25, 0.25, 0.25, 0.25], [0.0 1.0 2.0 1.0; 1.0 0.0 1.0 2.0; 2.0 1.0 0.0 1.0; 1.0 2.0 1.0 0.0]),
    ]

    for model in models
        Q = SeqSim.rate_matrix(model)
        π = SeqSim.get_frequencies(model)
        @test size(Q) == (4, 4)
        @test all(isapprox.(sum(Q, dims=1), 0.0; atol=1e-10))
        @test isapprox(Q * Diagonal(π), transpose(Q * Diagonal(π)); atol=1e-10)
        @test isapprox(Q * π, zeros(4); atol=1e-10)

        dec = SeqSim.decompose(Q)
        @test isapprox(Matrix(dec.V) * Diagonal(Vector(dec.λ)) * Matrix(dec.V⁻¹), Q; atol=1e-10)
    end
end

@testset "substitution model validation" begin
    @test_throws ArgumentError F81([0.1, 0.2, 0.3])
    @test_throws ArgumentError F81([0.1, 0.2, 0.3, 0.5])
    @test_throws ArgumentError F81([-0.1, 0.2, 0.3, 0.6])
    @test_throws ArgumentError K2P(0.0)
    @test_throws ArgumentError K2P(-1.0)
    @test_throws ArgumentError HKY([0.25, 0.25, 0.25, 0.25], 0.0)
    @test_throws ArgumentError GTR([0.25, 0.25, 0.25, 0.25], ones(3, 3))
    @test_throws ArgumentError GTR([0.25, 0.25, 0.25, 0.25], [0.0 1.0 2.0 1.0; 3.0 0.0 1.0 2.0; 2.0 1.0 0.0 1.0; 1.0 2.0 1.0 0.0])
    @test_throws ArgumentError SeqSim.rate_matrix([0.25, 0.25, 0.25, 0.25], [0.0 -1.0 1.0 1.0; -1.0 0.0 1.0 1.0; 1.0 1.0 0.0 1.0; 1.0 1.0 1.0 0.0])
end
