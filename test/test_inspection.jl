@testset "alignment inspection helpers" begin
    alignment = small_alignment_fixture()

    @test site_window(alignment, 2, 5) == [2, 3, 4, 5]
    @test centered_site_window(alignment, 1; radius=2) == [1, 2, 3]
    @test centered_site_window(alignment, 6; radius=2) == [4, 5, 6]
    @test centered_site_window(alignment, 3; radius=0) == [3]
    @test_throws ArgumentError site_window(alignment, 0, 2)
    @test_throws ArgumentError site_window(alignment, 4, 2)
    @test_throws ArgumentError centered_site_window(alignment, 7)
    @test_throws ArgumentError centered_site_window(alignment, 3; radius=-1)

    inspection = inspect_sites(alignment; sites=[6, 2], reference=1)
    @test inspection.sequence_count == 4
    @test inspection.site_count == 6
    @test inspection.selected_sites == [6, 2]
    @test inspection.reference == 1
    @test length(inspection.rows) == 2

    row6 = inspection.rows[1]
    @test row6.site == 6
    @test row6.observed_states == UInt8[1, 2, 4]
    @test row6.state_counts == [1, 2, 0, 1]
    @test row6.is_variable
    @test row6.minor_allele_count == 1
    @test row6.minor_allele_frequency == 0.25
    @test row6.is_singleton
    @test !row6.is_parsimony_informative
    @test row6.nonreference_frequency == 0.75

    row2 = inspection.rows[2]
    @test row2.site == 2
    @test row2.state_counts == [0, 2, 0, 2]
    @test row2.is_parsimony_informative
    @test row2.nonreference_frequency == 0.5
    @test [row.site for row in inspection.rows] == [6, 2]

    empty_inspection = inspect_sites(alignment; sites=Int[])
    @test empty_inspection.selected_sites == Int[]
    @test empty_inspection.rows == SiteInspectionRow[]

    window_inspection = inspect_window(alignment, 2, 3)
    @test window_inspection.selected_sites == [2, 3]
    centered = inspect_centered_window(alignment, 5; radius=1)
    @test centered.selected_sites == [4, 5, 6]

    projection = selected_site_strings(alignment; sites=[6, 2])
    @test projection == [
        (taxon = "a", value = "AC"),
        (taxon = "b", value = "CC"),
        (taxon = "c", value = "CT"),
        (taxon = "d", value = "TT"),
    ]
    @test selected_site_strings(alignment; sites=Int[]) == [
        (taxon = "a", value = ""),
        (taxon = "b", value = ""),
        (taxon = "c", value = ""),
        (taxon = "d", value = ""),
    ]

    compact = sprint(show, inspection)
    @test occursin("AlignmentInspection(4 sequences, 6 sites, 2 selected sites", compact)
    @test occursin("reference=1", compact)

    text = sprint(show, MIME"text/plain"(), inspection)
    @test occursin("AlignmentInspection", text)
    @test occursin("Alignment sites    : 6", text)
    @test occursin("Selected sites     : 2", text)
    @test occursin("Site coordinates   : original alignment positions", text)
    @test occursin("Selected-site preview (showing 2 of 2):", text)
    @test occursin("site 6", text)
    @test occursin("observed=ACT", text)
    @test !occursin("Source sites", text)
    @test !occursin("Site rows", text)
    @test !occursin("states=", text)
    @test occursin("nonref=0.75", text)

    full_inspection = inspect_sites(alignment; sites=:, reference=1)
    truncated = sprint(show, MIME"text/plain"(), full_inspection)
    @test occursin("Selected-site preview (showing 5 of 6):", truncated)
    @test occursin("... (1 more selected sites)", truncated)
    @test !occursin("more sites)", truncated)
end
