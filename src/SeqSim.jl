module SeqSim

using LinearAlgebra
using Parameters
using Phylo
using StatsBase
using UnPack

include("utils.jl")
include("substitution.jl")
include("sequencesimulator.jl")


export JC, F81, HKY, StrictClock, SiteModel
export nucleotides
export simulate_sequence, simulate_sequences!

end # module SeqSim
