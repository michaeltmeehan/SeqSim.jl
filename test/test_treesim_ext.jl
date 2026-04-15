push!(LOAD_PATH, normpath(@__DIR__, "..", "..", "EpiSim.jl"))
push!(LOAD_PATH, normpath(@__DIR__, "..", "..", "TreeSim.jl"))
using EpiSim
using TreeSim

function bridge_binary_tree()
    return Tree(
        [0.0, 0.4, 0.8, 1.1, 1.3],
        [2, 4, 0, 0, 0],
        [3, 5, 0, 0, 0],
        [0, 1, 1, 2, 2],
        [Root, Binary, SampledLeaf, SampledLeaf, SampledLeaf],
        [0, 0, 0, 0, 0],
        [900, 0, 103, 104, 105],
    )
end

function zero_rate_site_model(sequence_length::Int)
    return SiteModel(
        sequence_length,
        0.0,
        0,
        0.0,
        0.0,
        JC(),
        [collect(1:sequence_length)],
        [0.0],
    )
end

function fully_invariant_site_model(sequence_length::Int)
    return SiteModel(
        sequence_length,
        0.1,
        0,
        0.0,
        1.0,
        JC(),
        Vector{Vector{Int64}}(),
        Float64[],
    )
end

function composed_bridge_event_log()
    return EventLog(
        [0.0, 1.0, 2.0, 3.0, 4.0],
        [1, 2, 3, 2, 3],
        [0, 1, 1, 0, 0],
        [EK_Seeding, EK_Transmission, EK_Transmission, EK_SerialSampling, EK_SerialSampling],
    )
end

@testset "TreeSim extension bridge" begin
    @testset "canonical tree simulation returns layered outputs" begin
        tree = bridge_binary_tree()
        site_model = SiteModel(24, 0.15, 2, 1.0, 0.0, JC())
        result = simulate_tree_sequences(MersenneTwister(42), tree, site_model)

        @test result.root == root(tree)
        @test result.tips == tips(tree)
        @test length(result.node_sequences) == nnodes(tree)
        @test length(result.tip_alignment) == nleaves(tree)
        @test result.tip_alignment == result.node_sequences[result.tips]
        @test all(length(seq.value) == site_model.sequence_length for seq in result.node_sequences)
        @test all(seq.value == result.node_sequences[seq_index].value for (seq, seq_index) in zip(result.tip_alignment, result.tips))
    end

    @testset "node sequences are indexed directly by TreeSim node id" begin
        tree = bridge_binary_tree()
        site_model = SiteModel(MersenneTwister(3), 10, 0.1, 0, 0.0, 0.0, JC())
        result = simulate_tree_sequences(MersenneTwister(4), tree, site_model)

        @test eachindex(result.node_sequences) == eachindex(tree)
        for node in eachindex(tree)
            @test result.node_sequences[node].time == tree.time[node]
            expected_taxon = tree.label[node] == 0 ? "node_$node" : tree.label[node]
            @test result.node_sequences[node].taxon == expected_taxon
        end
    end

    @testset "metadata follows TreeSim node ids, labels, and times" begin
        tree = bridge_binary_tree()
        site_model = SiteModel(MersenneTwister(7), 8, 0.1, 0, 0.0, 0.0, HKY([0.25, 0.25, 0.25, 0.25], 2.0))
        result = simulate_tree_sequences(MersenneTwister(9), tree, site_model)

        @test result.node_sequences[1].taxon == 900
        @test result.node_sequences[2].taxon == "node_2"
        @test result.node_sequences[3].taxon == 103
        @test [seq.taxon for seq in result.tip_alignment] == [104, 105, 103]
        @test [seq.time for seq in result.node_sequences] == tree.time
        @test [seq.time for seq in result.tip_alignment] == tree.time[result.tips]
    end

    @testset "incoming branch lengths and traversal drive propagation reproducibly" begin
        tree = bridge_binary_tree()
        site_model = SiteModel(MersenneTwister(11), 32, 0.4, 0, 0.0, 0.0, F81([0.7, 0.1, 0.1, 0.1]))

        result = simulate_tree_sequences(MersenneTwister(13), tree, site_model)
        repeat = simulate_tree_sequences(MersenneTwister(13), tree, site_model)

        @test [seq.value for seq in result.node_sequences] == [seq.value for seq in repeat.node_sequences]
        @test branch_length(tree, 2) == tree.time[2] - tree.time[1]
        @test branch_length(tree, 4) == tree.time[4] - tree.time[2]
        @test result.node_sequences[4].value != result.node_sequences[2].value ||
              result.node_sequences[5].value != result.node_sequences[2].value ||
              result.node_sequences[3].value != result.node_sequences[1].value
    end

    @testset "public bridge wrappers are thin and reproducible" begin
        tree = bridge_binary_tree()
        site_model = SiteModel(MersenneTwister(17), 12, 0.2, 0, 0.0, 0.0, JC())

        tree_result = simulate_tree_sequences(MersenneTwister(19), tree, site_model)
        alias_result = simulate_on_tree(MersenneTwister(19), tree, site_model)
        tip_alignment = simulate_alignment(MersenneTwister(19), tree, site_model)

        @test [seq.value for seq in alias_result.node_sequences] == [seq.value for seq in tree_result.node_sequences]
        @test alias_result.tip_alignment == tree_result.tip_alignment
        @test tip_alignment == tree_result.tip_alignment
    end

    @testset "public wrapper reproducibility and tip projection" begin
        tree = bridge_binary_tree()
        site_model = SiteModel(MersenneTwister(23), 20, 0.15, 2, 1.0, 0.2, HKY([0.25, 0.25, 0.25, 0.25], 2.0))

        result_a = simulate_tree_sequences(MersenneTwister(29), tree, site_model)
        result_b = simulate_tree_sequences(MersenneTwister(29), tree, site_model)
        alignment_a = simulate_alignment(MersenneTwister(29), tree, site_model)
        alignment_b = simulate_alignment(MersenneTwister(29), tree, site_model)

        @test [seq.value for seq in result_a.node_sequences] == [seq.value for seq in result_b.node_sequences]
        @test result_a.tip_alignment == result_b.tip_alignment
        @test alignment_a == alignment_b
        @test result_a.tip_alignment == alignment_a
        @test result_a.tip_alignment == result_a.node_sequences[tips(tree)]
    end

    @testset "zero mutation rate preserves the root sequence through public wrappers" begin
        tree = bridge_binary_tree()
        site_model = zero_rate_site_model(16)

        result = simulate_tree_sequences(MersenneTwister(31), tree, site_model)
        alignment = simulate_alignment(MersenneTwister(31), tree, site_model)
        root_value = result.node_sequences[result.root].value

        @test all(seq.value == root_value for seq in result.node_sequences)
        @test all(seq.value == root_value for seq in result.tip_alignment)
        @test alignment == result.tip_alignment
    end

    @testset "full-invariance regime preserves all sites through public wrappers" begin
        tree = bridge_binary_tree()
        site_model = fully_invariant_site_model(18)

        result = simulate_tree_sequences(MersenneTwister(37), tree, site_model)
        alignment = simulate_alignment(MersenneTwister(37), tree, site_model)
        root_value = result.node_sequences[result.root].value

        @test isempty(site_model.variable_sites)
        @test isempty(site_model.μ)
        @test all(seq.value == root_value for seq in result.node_sequences)
        @test alignment == result.tip_alignment
    end

    @testset "public alignment wrapper preserves tip order and metadata" begin
        tree = bridge_binary_tree()
        site_model = SiteModel(MersenneTwister(41), 10, 0.1, 0, 0.0, 0.0, JC())

        result = simulate_tree_sequences(MersenneTwister(43), tree, site_model)
        alignment = simulate_alignment(MersenneTwister(43), tree, site_model)
        tip_ids = tips(tree)

        @test result.tips == tip_ids
        @test [seq.taxon for seq in alignment] == [tree.label[node] == 0 ? "node_$node" : tree.label[node] for node in tip_ids]
        @test [seq.time for seq in alignment] == tree.time[tip_ids]
        @test alignment == result.node_sequences[tip_ids]
    end

    @testset "zero branch lengths are outside the canonical TreeSim bridge contract" begin
        zero_edge_tree = Tree(
            [0.0, 0.0, 0.7],
            [2, 0, 0],
            [3, 0, 0],
            [0, 1, 1],
            [Root, SampledLeaf, SampledLeaf],
            [0, 0, 0],
            [0, 101, 102],
        )
        site_model = SiteModel(8, 0.1, 0, 0.0, 0.0, JC())

        @test_throws ArgumentError simulate_tree_sequences(MersenneTwister(47), zero_edge_tree, site_model)
        @test_throws ArgumentError simulate_alignment(MersenneTwister(47), zero_edge_tree, site_model)
    end

    @testset "invalid rooted trees are rejected at the bridge boundary" begin
        empty = Tree()
        @test_throws ArgumentError simulate_tree_sequences(MersenneTwister(1), empty, SiteModel(4, 0.1, 0, 0.0, 0.0, JC()))

        multi_root = Tree(
            [0.0, 0.0],
            [0, 0],
            [0, 0],
            [0, 0],
            [Root, Root],
            [0, 0],
            [1, 2],
        )
        site_model = SiteModel(4, 0.1, 0, 0.0, 0.0, JC())
        @test_throws ArgumentError simulate_tree_sequences(MersenneTwister(2), multi_root, site_model)
    end
end

@testset "Composed EpiSim to TreeSim to SeqSim bridge" begin
    @testset "event log sampled ancestry simulates a tip alignment" begin
        log = composed_bridge_event_log()
        tree = tree_from_eventlog(log)
        site_model = SiteModel(MersenneTwister(101), 14, 0.2, 0, 0.0, 0.0, JC())

        result = simulate_tree_sequences(MersenneTwister(102), tree, site_model)
        alignment = simulate_alignment(MersenneTwister(102), tree, site_model)

        @test validate_event_log(log; throw=false)
        @test validate_tree(tree)
        @test validate_tree_against_eventlog(log, tree)
        @test length(alignment) == nleaves(tree)
        @test alignment == result.tip_alignment
        @test result.tip_alignment == result.node_sequences[tips(tree)]
        @test [seq.taxon for seq in alignment] == tree.label[tips(tree)]
        @test [seq.time for seq in alignment] == tree.time[tips(tree)]
    end

    @testset "composed workflow is reproducible at the alignment level" begin
        log = composed_bridge_event_log()
        tree = tree_from_eventlog(log)
        site_model = SiteModel(MersenneTwister(201), 16, 0.3, 2, 1.0, 0.1, HKY([0.25, 0.25, 0.25, 0.25], 2.0))

        alignment_a = simulate_alignment(MersenneTwister(202), tree, site_model)
        alignment_b = simulate_alignment(MersenneTwister(202), tree, site_model)
        tree_result = simulate_tree_sequences(MersenneTwister(202), tree, site_model)

        @test alignment_a == alignment_b
        @test alignment_a == tree_result.tip_alignment
    end

    @testset "empty sampled ancestry is rejected by SeqSim simulation" begin
        log = EventLog(
            [0.0, 1.0, 2.0],
            [1, 2, 2],
            [0, 1, 0],
            [EK_Seeding, EK_Transmission, EK_Removal],
        )
        tree = tree_from_eventlog(log)
        site_model = SiteModel(8, 0.1, 0, 0.0, 0.0, JC())

        @test isempty(tree)
        @test validate_tree_against_eventlog(log, tree)
        @test_throws ArgumentError simulate_tree_sequences(MersenneTwister(301), tree, site_model)
        @test_throws ArgumentError simulate_alignment(MersenneTwister(301), tree, site_model)
    end
end
