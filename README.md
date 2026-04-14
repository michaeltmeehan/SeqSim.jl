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
