# SeqSim.jl

`SeqSim.jl` is a standalone sequence and alignment simulation package. It
provides a small, explicit DNA sequence/alignment core rather than a broad
phylogenetic workflow layer.

## Installation

From this repository:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

To use the package from another local Julia environment, add it by path:

```julia
using Pkg
Pkg.add(path="/path/to/outbreak-recovery/repos/SeqSim.jl")
```

## Supported scope

The active core currently supports:

- unambiguous DNA sequences over `A`, `C`, `G`, and `T`
- `Sequence` objects with optional `taxon` and `time` metadata
- alignments represented as `Vector{Sequence}` with equal sequence lengths
- root sequence simulation from base frequencies or a `SiteModel`
- sequence propagation over elapsed time using `SequencePropagator`
- sequence and site filtering utilities over equal-length alignments
- SNP/variable-site detection, SNP-only extraction, and reference-based SNP counts
- substitution models `JC`, `F81`, `K2P`, `HKY`, and `GTR`
- site models with invariant sites and optional discrete gamma rate categories
- optional tree-driven simulation from `TreeSim.Tree` via a package extension
- write-only alignment export to FASTA, NEXUS, and PHYLIP

## Explicitly unsupported for now

The current standalone core does not support:

- RNA, amino acids, ambiguity codes, gaps, or missing sequence states
- read-side sequence/alignment IO
- broad tree package integration or orchestration workflows
- broad format conversion
- modelling features beyond the substitution/site models listed above

Tree-facing behaviour is kept behind narrow optional adapters rather than by
making the core package depend on upstream tree packages.

## Core concepts

`Sequence` is the external interaction-facing sequence type. It stores a DNA
string plus optional `taxon` and `time` metadata. A supported alignment is a
`Vector{Sequence}` whose sequences are all the same length.

Internally, simulation and processing helpers use the efficient DNA state order
`A => 1`, `C => 2`, `G => 3`, and `T => 4`. User-facing helpers return
`Sequence` values, site indices, counts, matrices, or small inspection objects
rather than replacing the alignment representation with display-oriented data.

`AlignmentSummary` is a compact semantic summary of one alignment. It reports
sequence count and site classifications over alignment columns: variable sites,
invariant sites, singleton variable sites, and parsimony-informative sites.

## Minimal simulation example

```julia
using Random
using SeqSim

rng = MersenneTwister(1)
substitution_model = HKY([0.25, 0.25, 0.25, 0.25], 2.0)
site_model = SiteModel(12, 0.1, 2, 1.0, 0.0, substitution_model)

root = rand_seq(rng, site_model; taxon="root", time=0.0)
propagator = SequencePropagator(site_model)
child_states = propagator(rng, SeqSim.encode(root.value), 0.5)
child = Sequence(child_states; taxon="child", time=0.5)

alignment = [root, child]
write_alignment("alignment.fasta", alignment)
```

This creates a 12-site root sequence from the `SiteModel`, propagates it for
elapsed time `0.5`, stores the root and child as a two-sequence alignment, and
writes the alignment to FASTA. `rand_seq` returns a `Sequence`. The
`SequencePropagator` call returns encoded DNA states, so the example wraps them
back into a `Sequence` for external use.

The `SiteModel(12, 0.1, 2, 1.0, 0.0, substitution_model)` arguments are:
sequence length, mutation rate, number of discrete gamma categories, gamma
shape, invariant-site proportion, and substitution model. In this small example
the substitution model is `HKY` with equal base frequencies and
transition/transversion ratio `2.0`.

## Optional TreeSim bridge

When `TreeSim.jl` is loaded, SeqSim exposes extension methods for explicit
bridge functions:

```julia
using Random
using SeqSim
using TreeSim

tree = Tree(
    [0.0, 0.5, 0.7],
    [2, 0, 0],
    [3, 0, 0],
    [0, 1, 1],
    [Root, SampledLeaf, SampledLeaf],
    [0, 0, 0],
    [0, 101, 102],
)

site_model = SiteModel(12, 0.1, 0, 0.0, 0.0, JC())
simulation = simulate_tree_sequences(MersenneTwister(1), tree, site_model)
alignment = simulate_alignment(MersenneTwister(1), tree, site_model)

simulation.node_sequences # primary all-node result, indexed by TreeSim node id
simulation.tip_alignment  # derived tip-only Vector{Sequence}, ordered like simulation.tips
simulation.root           # TreeSim root id used for propagation
simulation.tips           # TreeSim tip ids in the same order as tip_alignment
alignment                 # the same tip-only product, returned directly
```

`simulate_alignment(rng, tree, site_model)` is a thin wrapper that returns only
the derived `tip_alignment`; use `simulate_tree_sequences` when all node
sequences are needed. `simulate_on_tree` remains available as a compatibility
alias for `simulate_tree_sequences`.

The bridge assumes a canonical rooted `TreeSim.Tree`: one root, reachable nodes,
parent-child consistency, preorder traversal, and TreeSim's node-id-as-vector-
index contract. SeqSim uses `TreeSim.branch_length(tree, node)` as the
nonnegative evolutionary branch length on the incoming edge to each non-root
node. The rooted structure and branch lengths drive simulation. TreeSim labels
and node times are preserved as sequence metadata: nonzero `tree.label[node]`
becomes `Sequence.taxon`, label `0` becomes `"node_$(node)"`, and
`tree.time[node]` becomes `Sequence.time`. Absolute node times are metadata in
the bridge result, not the direct evolutionary simulation primitive. TreeSim
remains an optional weak dependency.

### Composed EpiSim to TreeSim to SeqSim path

When both optional bridges are loaded, a sampled `EpiSim.EventLog` can be
converted to a TreeSim sampled-ancestry tree and then used for sequence
simulation:

```julia
using EpiSim
using Random
using SeqSim
using TreeSim

log = EventLog(
    [0.0, 1.0, 2.0, 3.0, 4.0],
    [1, 2, 3, 2, 3],
    [0, 1, 1, 0, 0],
    [EK_Seeding, EK_Transmission, EK_Transmission, EK_SerialSampling, EK_SerialSampling],
)

tree = tree_from_eventlog(log)
site_model = SiteModel(12, 0.1, 0, 0.0, 0.0, JC())

alignment = simulate_alignment(MersenneTwister(2), tree, site_model)
```

`tree_from_eventlog` extracts retained sampled ancestry, not the full outbreak
history. The resulting alignment is ordered like `TreeSim.tips(tree)`, with tip
labels and times preserved through the SeqSim metadata policy.

## Small alignment example

The repository includes the same example as `examples/small_alignment.jl`. It is
small enough to inspect by eye and is reused as the test fixture for processing,
inspection, and lightweight viewer checks.

```julia
using SeqSim

alignment = [
    Sequence("ACGTAA"; taxon="a"),
    Sequence("ACGTAC"; taxon="b"),
    Sequence("ATGTCC"; taxon="c"),
    Sequence("ATGTCT"; taxon="d"),
]

AlignmentSummary(alignment)
```

The alignment has four sequences and six sites. Read each site as a vertical
column. Sites 1, 3, and 4 are invariant because every sequence has the same
base. Sites 2, 5, and 6 are variable because at least two bases are observed.
Sites 2 and 5 are parsimony-informative because two states each appear in at
least two sequences. Site 6 is a singleton variable site because at least one
observed state appears only once.

`AlignmentSummary(alignment)` counts those classes across all sites in the
alignment. It does not create a new alignment representation.

## Processing and SNP utilities

Processing helpers operate on the external `Alignment` view
(`Vector{Sequence}`) and preserve the simulation-facing representation. They
return filtered `Sequence` copies or plain site/count vectors, rather than
introducing display-oriented storage.

```julia
variable_sites(alignment)       # [2, 5, 6]
invariant_sites(alignment)      # [1, 3, 4]
snp_alignment(alignment)        # new alignment restricted to sites 2, 5, and 6
filter_sites(alignment, [6, 2]) # preserves requested site order
```

`variable_sites` and `invariant_sites` return original alignment site
coordinates. `snp_alignment` returns a filtered `Vector{Sequence}` containing
only the variable sites. `filter_sites` accepts `:`, integer vectors, integer
ranges, or boolean masks; integer vectors preserve caller-specified order, so
`[6, 2]` projects site 6 first and site 2 second.

Reference-based SNP helpers count mismatches from a reference sequence. The
reference can be a sequence index, a `Sequence`, or a DNA string.

```julia
site_snp_counts(alignment; reference=1)
# [0, 2, 0, 0, 2, 3]

sequence_snp_counts(alignment; reference=1)
# [0, 1, 3, 3]
```

With `reference=1`, the first sequence is the reference. `site_snp_counts`
returns one count per alignment site: at site 6, three of the four sequences
differ from sequence 1. `sequence_snp_counts` returns one count per sequence:
sequence 1 has zero differences from itself, while sequences 3 and 4 differ at
three sites.

Site-count helpers use the computation-facing DNA state order `A, C, G, T`:

```julia
site_state_counts(alignment; sites=[2, 6])
# 4 x 2 Matrix:
# 0  1
# 2  2
# 0  0
# 2  1
```

The first selected site is site 2, where two sequences have `C` and two have
`T`. The second selected site is site 6, where the counts are one `A`, two `C`,
and one `T`.

Other useful derived outputs:

```julia
minor_allele_frequencies(alignment) # [0.0, 0.5, 0.0, 0.0, 0.5, 0.25]
pairwise_difference_matrix(alignment)
consensus_sequence(alignment)       # Sequence("ACGTAC"; taxon="consensus")
```

Because the current core only supports unambiguous DNA without gaps or
ambiguity codes, SNP sites and variable sites are equivalent for now. Gap and
ambiguous-site helpers remain deferred until those alphabets are explicitly
supported.

Site selectors are intentionally narrow: use `:`, integer vectors, integer
ranges, or boolean masks. Duplicate integer sites are rejected; caller-specified
order is preserved for integer selectors. Empty selectors are supported by
derived numeric/inspection helpers where zero-site outputs are meaningful, but
`filter_sites` rejects them because `Sequence` values must be non-empty.
Reference-based count helpers accept a sequence index, a `Sequence`, or a DNA
string. By default, both `site_snp_counts` and `sequence_snp_counts` count
across all sites; pass `sites=variable_sites(alignment)` for variable-site-only
counts.

## Inspection helpers

Inspection helpers are derived view-models for scripts and future viewer work.
They preserve original alignment site coordinates and do not replace alignment
storage.

```julia
site_window(alignment, 2, 5)                  # [2, 3, 4, 5]
centered_site_window(alignment, 1; radius=2)  # [1, 2, 3]

inspection = inspect_sites(alignment; sites=variable_sites(alignment), reference=1)
inspection.rows[1].site
inspection.rows[1].minor_allele_frequency
inspection.rows[1].nonreference_frequency

inspect_window(alignment, 2, 4)
selected_site_strings(alignment; sites=[6, 2])
# [
#   (taxon = "a", value = "AC"),
#   (taxon = "b", value = "CC"),
#   (taxon = "c", value = "CT"),
#   (taxon = "d", value = "TT"),
# ]
```

`inspect_sites` returns an `AlignmentInspection`: one row per selected original
site, with observed encoded states, counts in `A,C,G,T` order, variable-site
flags, minor-allele frequency, singleton/parsimony-informative flags, and
optional nonreference frequency. `selected_sites` and `row.site` use original
alignment positions; sites are not reindexed after selection. The `reference`
argument is used only for reference-aware values such as nonreference
frequency.

## Lightweight viewer

The lightweight viewer is a bounded textual layer built on `AlignmentInspection`.
It is intended for REPL, notebook, and script previews; it is not a GUI, plotting
system, or alternate alignment representation.

```julia
viewer = AlignmentViewer(alignment; sites=variable_sites(alignment),
                         reference=1, max_sequences=4, max_sites=6)
print(viewer_preview(viewer))

window_viewer(alignment, 2, 5)
centered_window_viewer(alignment, 5; radius=1)
variable_site_viewer(alignment)
```

Viewer previews preserve original site coordinates in their site headers,
truncate displayed sequences/sites to the configured limits, and can include a
small reference-aware row with mismatch markers. They intentionally avoid ANSI
color, plotting dependencies, GUI state, and tree-facing behavior.

## Validation

The current hardening target is described in `VALIDATION.md` and
`TRUST_CRITERIA.md`. Tests focus on representation invariants, fixed-seed
reproducibility, compact statistical checks, propagation invariants, and
write-side IO integrity for the supported formats.

## Testing

From the package directory:

```julia
using Pkg
Pkg.test()
```
