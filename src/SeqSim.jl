"""
Module `SeqSim` for simulating DNA sequences along the branches of a phylogenetic tree.

# Exports
- `JC, F81, K2P, HKY, GTR, rate_matrix`
- `StrictClock`
- `SiteModel, assign_rate_categories`
- `nucleotides, weights`
- `decompose, mod_wrap`
- `simulate_sequence, simulate_sequences!, compute_transition_weights!, propagate_sequence, assign_rates, update_site`
- `tip_sequences`
"""
module SeqSim

using BioSequences
using BioSymbols
using Distributions
using Lazy
using LinearAlgebra
using Parameters
using Phylo
using StatsBase
using TreeSim
using UnPack

include("utils.jl")
include("substitution.jl")
include("clock.jl")
include("site.jl")
include("sequencesimulator.jl")
include("sequenceanalysis.jl")


export JC, F81, K2P, HKY, GTR, rate_matrix 
export StrictClock
export SiteModel, assign_rate_categories
export nucleotides, weights
export decompose, mod_wrap
export simulate_sequences!, compute_transition_weights!, propagate_sequence, assign_rates, update_site
export tip_sequences

end # module SeqSim
