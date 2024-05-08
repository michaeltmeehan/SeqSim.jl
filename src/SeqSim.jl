module SeqSim

using BioSequences
using Distributions
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


export JC, F81, K2P, HKY, GTR, StrictClock
export SiteModel, site_rates, discretize_dist, rate_cat
export nucleotides, weights
export decompose, rate_matrix, simulate_sequence, simulate_sequences!, mod_wrap, compute_transition_weights!, propagate_sequence
export tip_sequences

end # module SeqSim
