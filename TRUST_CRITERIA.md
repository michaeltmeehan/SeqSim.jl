# Trust Criteria

## Role of this package

`SeqSim.jl` is the standalone sequence and alignment simulation package in the recovered outbreak-modelling ecosystem.

Its current role is to provide a trustworthy minimal sequence/alignment simulation core within a clearly documented scope.

---

## Trust goal

`SeqSim.jl` is trustworthy when its current simulation and representation core can be relied upon to generate and preserve sequences/alignments correctly within the documented supported domain.

---

## What trust does mean here

Trust in this phase means:

- the supported biological and technical scope is explicit
- the current sequence/alignment representations are stable
- simple simulation behaviour agrees with documented model expectations
- fixed-seed reproducibility works as documented
- supported IO paths preserve content correctly
- invalid or unsupported scenarios fail clearly

---

## What trust does not mean here

Trust in this phase does **not** imply:

- support for all phylogenetic workflows
- strong coupling to tree packages
- support for all substitution/site models
- support for all sequence formats
- orchestration-level workflow integration

---

## Stable trust boundary

The current trust boundary includes:

- the recovered minimal sequence/alignment simulation core
- sequence representation decisions currently retained in the active core
- alignment and IO behaviour intentionally exposed in the current minimal API

Anything outside this boundary should be treated as provisional unless explicitly promoted into the stable API.

---

## Conditions required for trust

### 1. Scope clarity

The package must clearly state what it simulates and what it does not.

### 2. Representation reliability

Core sequence/alignment objects must behave predictably and consistently.

### 3. Simulation reliability

Simple statistical checks must support the intended model behaviour.

### 4. Reproducibility

Fixed-seed simulation must behave deterministically where promised.

### 5. IO reliability

Supported IO must preserve data and metadata as documented.

---

## Known trust risks

Current or likely risks include:

- plausible-looking but quantitatively incorrect simulation output
- hidden issues in compact sequence representations
- metadata drift or corruption across IO paths
- unsupported biological assumptions being inferred by users
- insufficient edge-case coverage

---

## Required evidence before calling this package trustworthy

The following evidence is required:

- a written scope/model note
- representation invariant tests
- statistical validation tests for simple simulation settings
- fixed-seed reproducibility tests
- IO round-trip tests for supported formats
- explicit failure tests for invalid inputs

---

## Phase-2 completion standard

For the purposes of this project phase, `SeqSim.jl` is trustworthy when:

1. its supported scope is clearly stated
2. its core representations are stable and tested
3. its current simulation behaviour has basic quantitative validation
4. its IO behaviour is tested where supported
5. it is reliable enough to serve as the standalone sequence/alignment layer for the current ecosystem phase
