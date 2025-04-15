module SeqSim

using Crayons
using Distributions
using LinearAlgebra
# using LoopVectorization
using Random
using StaticArrays
using StatsBase


include("sequence.jl")
include("substitution.jl")
include("site.jl")
include("sequence_simulator.jl")


export JC, F81, K2P, HKY, GTR, SubstitutionModel
export SiteModel
export rand_seq, update_sequence!
export Sequence

end # module SeqSim
