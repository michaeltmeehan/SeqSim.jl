module SeqSim

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


export JC, F81, K2P, HKY, GTR, StrictClock, SiteModel
export nucleotides
export decompose, rate_matrix, simulate_sequence, simulate_sequences!, mod_wrap, compute_transition_weights!, propagate_sequence

end # module SeqSim
