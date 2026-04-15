@testset "lightweight alignment viewer" begin
    alignment = small_alignment_fixture()

    viewer = AlignmentViewer(alignment; max_sequences=2, max_sites=3)
    @test viewer.inspection.selected_sites == [1, 2, 3, 4, 5, 6]
    @test viewer.max_sequences == 2
    @test viewer.max_sites == 3
    @test !viewer.variable_only

    explicit = alignment_viewer(alignment; sites=[6, 2, 5], max_sequences=2, max_sites=2)
    @test explicit.inspection.selected_sites == [6, 2, 5]
    @test SeqSim.visible_sites(explicit) == [6, 2, 5]
    @test SeqSim.displayed_sites(explicit) == [6, 2]

    windowed = window_viewer(alignment, 2, 4)
    @test windowed.inspection.selected_sites == [2, 3, 4]

    centered = centered_window_viewer(alignment, 1; radius=2)
    @test centered.inspection.selected_sites == [1, 2, 3]

    variable = variable_site_viewer(alignment; sites=[6, 4, 2, 1])
    @test variable.inspection.selected_sites == [6, 2]
    @test SeqSim.visible_sites(variable) == [6, 2]

    variable_option = AlignmentViewer(alignment; sites=[1, 2, 5, 6], variable_only=true)
    @test SeqSim.visible_sites(variable_option) == [2, 5, 6]

    preview = viewer_preview(explicit)
    @test occursin("AlignmentViewer", preview)
    @test occursin("Sequences      : 4 (showing 2)", preview)
    @test occursin("Selected sites : 3 (showing 2)", preview)
    @test occursin("Sites: 6 2", preview)
    @test occursin("a           AC", preview)
    @test occursin("b           CC", preview)
    @test occursin("... (2 more sequences)", preview)
    @test occursin("... (1 more selected sites)", preview)

    ref_viewer = AlignmentViewer(alignment; sites=[6, 2], reference=1, max_sequences=2)
    ref_preview = viewer_preview(ref_viewer)
    @test occursin("Reference      : 1", ref_preview)
    @test occursin("Nonref: 0.75 0.5", ref_preview)
    @test occursin("reference   AC", ref_preview)
    @test occursin("            ..", ref_preview)
    @test occursin("            ^.", ref_preview)

    compact = sprint(show, ref_viewer)
    @test occursin("AlignmentViewer(4 sequences, 2 selected sites", compact)
    @test occursin("reference=1", compact)

    text = sprint(show, MIME"text/plain"(), variable_option)
    @test occursin("Site filter    : variable only", text)
    @test occursin("Sites: 2 5 6", text)

    empty = AlignmentViewer(alignment; sites=Int[])
    empty_preview = viewer_preview(empty)
    @test occursin("Selected sites : 0 (showing 0)", empty_preview)
    @test occursin("Sites          : <none>", empty_preview)

    @test_throws ArgumentError AlignmentViewer(alignment; max_sequences=0)
    @test_throws ArgumentError AlignmentViewer(alignment; max_sites=0)
end
