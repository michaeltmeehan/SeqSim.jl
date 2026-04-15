@testset "processing and SNP utilities" begin
    alignment = [
        Sequence("ACGTAC"; taxon="ref", time=0.0),
        Sequence("ATGTTC"; taxon="sample-1", time=1.0),
        Sequence("ACCTAA"; taxon="sample-2", time=2.0),
    ]

    @test variable_sites(alignment) == [2, 3, 5, 6]
    @test get_snps(alignment) == [2, 3, 5, 6]
    @test invariant_sites(alignment) == [1, 4]
    @test site_states(alignment, 5) == UInt8[1, 4, 1]

    site_filtered = filter_sites(alignment, [6, 2])
    @test [seq.value for seq in site_filtered] == ["CC", "CT", "AC"]
    @test site_filtered[1].taxon == "ref"
    @test site_filtered[2].time == 1.0
    @test [seq.value for seq in filter_sites(alignment, 2:4)] == ["CGT", "TGT", "CCT"]
    @test [seq.value for seq in filter_sites(alignment, [false, true, false, false, true, false])] == ["CA", "TT", "CA"]
    @test [seq.value for seq in filter_sites(alignment, :)] == [seq.value for seq in alignment]

    sequence_filtered = filter_sequences(seq -> seq.time !== nothing && seq.time >= 1.0, alignment)
    @test [seq.taxon for seq in sequence_filtered] == ["sample-1", "sample-2"]

    snps = snp_alignment(alignment)
    @test [seq.value for seq in snps] == ["CGAC", "TGTC", "CCAA"]
    @test invariant_filtered_alignment(alignment) == snps

    @test site_snp_counts(alignment) == [0, 1, 1, 0, 1, 1]
    @test site_snp_counts(alignment; sites=Int[]) == Int[]
    @test site_snp_counts(alignment; sites=[false, false, false, false, false, false]) == Int[]
    @test site_snp_counts(alignment; sites=[2, 5]) == [1, 1]
    @test site_snp_counts(alignment; reference=2) == [0, 2, 1, 0, 2, 1]
    @test site_snp_counts(alignment; reference=Sequence("ACGTTC")) == [0, 1, 1, 0, 2, 1]
    @test site_snp_counts(alignment; reference="ACGTTC") == [0, 1, 1, 0, 2, 1]
    @test sequence_snp_counts(alignment) == [0, 2, 2]
    @test sequence_snp_counts(alignment; sites=variable_sites(alignment)) == [0, 2, 2]
    @test sequence_snp_counts(alignment; reference="ATGTTC") == [2, 0, 4]
    @test sequence_snp_counts(alignment; sites=Int[]) == [0, 0, 0]
    @test sequence_snp_counts(alignment; sites=[2, 5]) == [0, 2, 0]
    @test sequence_snp_counts(alignment; reference=2, sites=:) == [2, 0, 4]

    invariant_alignment = [Sequence("ACGT"; taxon="a"), Sequence("ACGT"; taxon="b")]
    @test variable_sites(invariant_alignment) == Int[]
    @test invariant_sites(invariant_alignment) == [1, 2, 3, 4]
    @test sequence_snp_counts(invariant_alignment) == [0, 0]
    @test_throws ArgumentError snp_alignment(invariant_alignment)

    @test_throws ArgumentError filter_sites(alignment, Int[])
    @test_throws ArgumentError filter_sites(alignment, [0])
    @test_throws ArgumentError filter_sites(alignment, [7])
    @test_throws ArgumentError filter_sites(alignment, [1, 1])
    @test_throws ArgumentError filter_sites(alignment, [true, false])
    @test_throws ArgumentError filter_sites(alignment, [false, false, false, false, false, false])
    @test_throws ArgumentError site_states(alignment, 7)
    @test_throws ArgumentError site_snp_counts(alignment; reference=0)
    @test_throws ArgumentError site_snp_counts(alignment; reference=Sequence("ACG"))
    @test_throws ArgumentError site_snp_counts(alignment; reference="ACGN")
    @test_throws ArgumentError site_snp_counts(alignment; reference="ACG")
    @test_throws ArgumentError sequence_snp_counts(alignment; reference=4)
end

@testset "site counts, classifications, differences, and summaries" begin
    alignment = [
        Sequence("ACGTAA"; taxon="a"),
        Sequence("ACGTAC"; taxon="b"),
        Sequence("ATGTCC"; taxon="c"),
        Sequence("ATGTCT"; taxon="d"),
    ]

    counts = site_state_counts(alignment)
    @test size(counts) == (4, 6)
    @test counts[:, 1] == [4, 0, 0, 0]
    @test counts[:, 2] == [0, 2, 0, 2]
    @test counts[:, 5] == [2, 2, 0, 0]
    @test counts[:, 6] == [1, 2, 0, 1]
    @test site_state_counts(alignment; sites=[6, 2])[:, 1] == [1, 2, 0, 1]
    @test size(site_state_counts(alignment; sites=Int[])) == (4, 0)

    observed = observed_state_counts(alignment; sites=[2, 6])
    @test observed == [[0, 2, 0, 2], [1, 2, 0, 1]]

    @test variable_sites(alignment) == [2, 5, 6]
    @test invariant_sites(alignment) == [1, 3, 4]
    @test minor_allele_counts(alignment) == [0, 2, 0, 0, 2, 1]
    @test minor_allele_frequencies(alignment) == [0.0, 0.5, 0.0, 0.0, 0.5, 0.25]
    @test nonreference_frequencies(alignment; reference=1) == [0.0, 0.5, 0.0, 0.0, 0.5, 0.75]
    @test singleton_sites(alignment) == [6]
    @test parsimony_informative_sites(alignment) == [2, 5]

    @test pairwise_differences(alignment[1], alignment[2]) == 1
    @test pairwise_differences(alignment[1], alignment[2]; sites=Int[]) == 0
    @test pairwise_differences(alignment[1], alignment[4]; sites=[2, 5, 6]) == 3
    @test_throws ArgumentError pairwise_differences(Sequence("ACGT"), Sequence("ACG"))

    matrix = pairwise_difference_matrix(alignment)
    @test matrix == [0 1 3 3; 1 0 2 3; 3 2 0 1; 3 3 1 0]
    @test pairwise_difference_matrix(alignment; sites=[2, 5]) == [0 0 2 2; 0 0 2 2; 2 2 0 0; 2 2 0 0]

    consensus = consensus_sequence(alignment)
    @test consensus.value == "ACGTAC"
    @test consensus.taxon == "consensus"
    @test_throws ArgumentError consensus_sequence(alignment; tie=:error)
    @test_throws ArgumentError consensus_sequence(alignment; tie=:ambiguous)

    summary = AlignmentSummary(alignment)
    @test summary.sequence_count == 4
    @test summary.site_count == 6
    @test summary.variable_site_count == 3
    @test summary.invariant_site_count == 3
    @test summary.singleton_site_count == 1
    @test summary.parsimony_informative_site_count == 2
end
