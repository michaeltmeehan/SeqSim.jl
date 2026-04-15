# SeqSim.jl

`SeqSim.jl` is the standalone sequence and alignment simulation package in the
recovered outbreak-modelling ecosystem. In the current phase it provides a
small, explicit DNA sequence/alignment core rather than a broad phylogenetic
workflow layer.

## Supported scope

The active core currently supports:

- unambiguous DNA sequences over `A`, `C`, `G`, and `T`
- `Sequence` objects with optional `taxon` and `time` metadata
- alignments represented as `Vector{Sequence}` with equal sequence lengths
- root sequence simulation from base frequencies or a `SiteModel`
- sequence propagation over elapsed time using `SequencePropagator`
- sequence and site filtering utilities over validated alignments
- SNP/variable-site detection, SNP-only extraction, and reference-based SNP counts
- substitution models `JC`, `F81`, `K2P`, `HKY`, and `GTR`
- site models with invariant sites and optional discrete gamma rate categories
- write-only alignment export to FASTA, NEXUS, and PHYLIP

## Explicitly unsupported for now

The current standalone core does not support:

- RNA, amino acids, ambiguity codes, gaps, or missing sequence states
- read-side sequence/alignment IO
- tree package integration or orchestration workflows
- broad format conversion
- modelling features beyond the substitution/site models listed above

Future tree-facing behaviour should be added through narrow adapters rather
than by making this package depend on upstream tree packages.

## Basic usage

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

## Processing and SNP utilities

Processing helpers operate on the external `Alignment` view
(`Vector{Sequence}`) and preserve the simulation-facing representation. They
return filtered `Sequence` copies or plain site/count vectors, rather than
introducing display-oriented storage.

```julia
alignment = [
    Sequence("ACGTAC"; taxon="ref"),
    Sequence("ATGTTC"; taxon="sample-1"),
    Sequence("ACCTAA"; taxon="sample-2"),
]

variable_sites(alignment)          # [2, 3, 5, 6]
invariant_sites(alignment)         # [1, 4]
snp_alignment(alignment)           # alignment restricted to variable sites
site_snp_counts(alignment)         # non-reference counts by site
sequence_snp_counts(alignment)     # SNP counts by sequence
filter_sites(alignment, [6, 2])    # preserves requested site order
site_state_counts(alignment)       # 4 x site count matrix in A,C,G,T order
minor_allele_frequencies(alignment)
pairwise_difference_matrix(alignment)
consensus_sequence(alignment)
AlignmentSummary(alignment)
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

The computation-facing DNA state order is `A => 1`, `C => 2`, `G => 3`, and
`T => 4`. `site_states` and `site_state_counts` use this encoded order rather
than display-oriented characters.

## Inspection helpers

Inspection helpers are derived view-models for scripts and future viewer work.
They preserve original alignment site coordinates and do not replace alignment
storage.

```julia
site_window(alignment, 10, 25)
centered_site_window(alignment, 42; radius=5)

inspection = inspect_sites(alignment; sites=variable_sites(alignment), reference=1)
inspection.rows[1].site                 # original alignment site index
inspection.rows[1].minor_allele_frequency
inspection.rows[1].nonreference_frequency

inspect_window(alignment, 10, 25)
selected_site_strings(alignment; sites=[12, 4, 8])
```

## Lightweight viewer

The lightweight viewer is a bounded textual layer built on `AlignmentInspection`.
It is intended for REPL, notebook, and script previews; it is not a GUI, plotting
system, or alternate alignment representation.

```julia
viewer = AlignmentViewer(alignment; sites=variable_sites(alignment),
                         reference=1, max_sequences=10, max_sites=40)
viewer_preview(viewer)

window_viewer(alignment, 10, 25)
centered_window_viewer(alignment, 42; radius=5)
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
