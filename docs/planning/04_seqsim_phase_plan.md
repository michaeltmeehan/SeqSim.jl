## `SeqSim.jl` next-phase plan

### Role

`SeqSim.jl` is the standalone sequence and alignment simulation core.

### Current position

`SeqSim.jl` has been recovered and validated within its current standalone simulation scope.

### Immediate priorities

#### 1. Processing utilities

- sequence filtering
- alignment filtering
- site filtering
- gap filtering
- ambiguous-site handling helpers
- conversion utilities between internal and user-facing forms

#### 2. SNP analysis / filtering

- identify variable sites
- identify SNP-only sites
- filter invariant sites
- extract SNP matrices
- summarize SNP counts by sequence or site
- site-level summary reports

#### 3. Alignment viewer

- readable alignment display
- compact summaries for large alignments
- highlighting of variable sites / SNPs
- optional annotation display
- site-range slicing for inspection

#### 4. Visualization support

- site variability plots
- SNP density summaries
- sequence divergence summaries
- composition plots

### Medium-term priorities

- summary objects
- ensemble support
- export / conversion helpers

### Explicit deferrals

- direct tree-native identity
- upstream dependence on tree structures
- orchestration-oriented workflow logic
- broad fitting / calibration layers

### Architectural guardrails

- preserve efficient internal sequence and alignment encodings
- distinguish internal computation-facing types from external interaction-facing views
- build viewers and summaries on top of the internal core
- favor explicit conversion utilities rather than hidden coercions
