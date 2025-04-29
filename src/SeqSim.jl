module SeqSim

using Crayons
using Distributions
using LinearAlgebra
# using LoopVectorization
using Random
using StaticArrays
using StatsBase


include("Sequence.jl")
include("substitution.jl")
include("site.jl")
include("sequence_simulator.jl")
include("SequencePropagator.jl")
include("export.jl")


export JC, F81, K2P, HKY, GTR, SubstitutionModel
export SiteModel
export rand_seq, update_sequence!, SequencePropagator
export Sequence
export write_fasta, write_nexus, write_phylip, write_alignment

end # module SeqSim
