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
    @test_throws ArgumentError simulate_sequence(-1, [0.25, 0.25, 0.25, 0.25])
    @test_throws ArgumentError simulate_sequence(10, [0.25, 0.25, 0.25])
    @test_throws ArgumentError simulate_sequence(10, [0.25, 0.25, 0.25, -0.25])
    @test_throws ArgumentError simulate_sequence(10, [0.3, 0.3, 0.3, 0.3])
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

