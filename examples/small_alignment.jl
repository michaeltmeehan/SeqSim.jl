using SeqSim

# A tiny alignment for inspection and viewer examples.
#
# Read columns vertically:
# - sites 1, 3, and 4 are invariant
# - sites 2 and 5 are parsimony-informative variable sites
# - site 6 is a singleton variable site
alignment = [
    Sequence("ACGTAA"; taxon="a"),
    Sequence("ACGTAC"; taxon="b"),
    Sequence("ATGTCC"; taxon="c"),
    Sequence("ATGTCT"; taxon="d"),
]

println("Alignment")
for seq in alignment
    println(rpad(string(seq.taxon), 8), seq.value)
end
println()

# AlignmentSummary counts site classes across the whole alignment.
summary = AlignmentSummary(alignment)
println("Whole-alignment summary")
show(stdout, MIME"text/plain"(), summary)
println()
println()

# Site classification helpers return original alignment coordinates.
variables = variable_sites(alignment)
invariants = invariant_sites(alignment)
singletons = singleton_sites(alignment)
informative = parsimony_informative_sites(alignment)

println("Site classes")
println("  variable sites              : ", variables)
println("  invariant sites             : ", invariants)
println("  singleton variable sites    : ", singletons)
println("  parsimony-informative sites : ", informative)
println()

# Inspection rows are a derived view over selected original sites.
# The reference argument is used only for nonreference frequencies.
inspection = inspect_sites(alignment; sites=variables, reference=1)
println("Inspection of selected variable sites")
show(stdout, MIME"text/plain"(), inspection)
println()
println()

# Projections preserve the requested selected-site order.
println("Selected-site strings for sites [6, 2]")
for projection in selected_site_strings(alignment; sites=[6, 2])
    println(rpad(string(projection.taxon), 8), projection.value)
end
println()

# The lightweight viewer is a bounded text preview, not a new alignment type.
viewer = variable_site_viewer(alignment; reference=1)
println("Lightweight viewer preview")
print(viewer_preview(viewer))
