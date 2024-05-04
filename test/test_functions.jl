@testset "simulate_sequence tests" begin
    # Test basic functionality
    seq = simulate_sequence(10, [0.25, 0.25, 0.25, 0.25])
    @test length(seq) == 10
    @test all(n -> n in nucleotides, seq)

    # Test with non-uniform distribution
    seq = simulate_sequence(10_000, [0.1, 0.2, 0.3, 0.4])
    counts = countmap(seq)
    @test counts['a'] < counts['t']
    @test counts['c'] < counts['g']

    # Test handling of incorrect inputs
    @test_throws ArgumentError simulate_sequence(10, [0.1, 0.2, 0.3])  # Not enough frequencies
    @test_throws ArgumentError simulate_sequence(10, [0.1, 0.2, 0.3, 0.5])  # Sum != 1
end