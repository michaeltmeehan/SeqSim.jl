# Validation Plan

## Purpose

This document defines the validation work required for `SeqSim.jl` to be considered a trustworthy standalone package for sequence and alignment simulation.

`SeqSim.jl` is intended to remain standalone for now. Any future tree-facing behaviour should occur through narrow adapters rather than redefining the package around upstream dependencies.

---

## Current scope

`SeqSim.jl` currently supports:

- a recovered and stabilised sequence/alignment simulation core
- clarified sequence representation boundaries
- a trimmed dependency surface
- a minimal loadable standalone package structure

---

## Out of scope

The following are explicitly out of scope for the current phase:

- strong coupling to upstream tree packages
- orchestration-layer workflows
- full phylogenetic workflow integration
- unsupported substitution/site-model features not explicitly documented
- broad format/conversion ecosystems beyond the current active core

---

## Core validation questions

Validation work must establish that:

1. the current sequence simulation core behaves in line with its documented models
2. sequence and alignment representations are stable and semantically clear
3. seeded simulation is reproducible
4. IO, where supported, preserves content and metadata correctly
5. unsupported biological or technical scenarios are clearly excluded

---

## Required validation areas

### 1. Scope and model specification

Validation must begin by documenting:

- supported sequence/alphabet types
- supported alignment assumptions
- supported substitution/site models
- supported metadata assumptions
- explicitly unsupported features

### 2. Statistical validation of simulation behaviour

Tests must examine:

- empirical base/state frequencies where applicable
- expected transition behaviour for simple models
- sequence-length invariants
- site-level behaviour in simple controlled settings
- reproducibility under fixed seeds

### 3. Representation validation

Tests must verify:

- stable handling of the core sequence representation
- correct length/content semantics
- correct identifier/metadata association where applicable
- rejection of invalid symbols or malformed inputs where appropriate

### 4. IO round-trip validation

Where IO is supported, tests must verify:

- write-read round-trip preservation of sequence content
- preservation of identifiers
- preservation of alignment dimensions
- preservation of metadata where intended
- failure behaviour for invalid or unsupported inputs

### 5. Regression and edge-case validation

Tests should include:

- tiny sequences
- tiny alignments
- unusual but valid lengths
- invalid symbols
- metadata mismatch cases
- deterministic fixed-seed regression examples

---

## Canonical benchmark suite

A benchmark suite should include:

- simple model simulations with expected frequency behaviour
- fixed-seed deterministic cases
- IO round-trip examples
- valid edge-case examples
- invalid examples expected to fail

Each benchmark should document:

- model/settings
- RNG setup
- input object(s)
- expected property or expected failure

---

## Downstream validation relevance

`SeqSim.jl` remains standalone in this phase.

Validation should still document what future adapters may safely assume about:

- sequence representation
- alignment representation
- metadata stability
- scope limitations

This should be stated as a stable local contract, not as an integration commitment.

---

## Exit criteria for Phase 2

`SeqSim.jl` should not be considered fully validated for this phase until:

- supported models and representations are explicitly documented
- simulation behaviour has basic statistical validation
- reproducibility tests exist
- IO round-trip tests exist for supported formats
- invalid/unsupported inputs fail clearly
- the minimal stable public API is clearly defined

---

## Evidence of successful validation

Evidence should include:

- statistical validation tests for simple models
- representation invariant tests
- IO round-trip tests
- fixed-seed reproducibility tests
- documented expected failures for unsupported inputs
