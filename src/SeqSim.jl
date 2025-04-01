module SeqSim

using Distributions
using LinearAlgebra
using Random
using StaticArrays
using StatsBase

include("substitution.jl")
include("site.jl")
include("sequence_simulator.jl")


export JC, F81, K2P, HKY, GTR, SubstitutionModel
export SiteModel
export update_sequence!

end # module SeqSim
