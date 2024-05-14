module SeqSim

using BioSequences
using Distributions
using Lazy
using LinearAlgebra
using Parameters
using Phylo
using StatsBase
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
export simulate_sequence, simulate_sequences!, compute_transition_weights!, propagate_sequence, assign_rates, update_site
export tip_sequences

end # module SeqSim
