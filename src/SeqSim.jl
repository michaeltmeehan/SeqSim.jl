module SeqSim

using Crayons
using Distributions
using LinearAlgebra
using Random
using StaticArrays


include("sequence.jl")
include("substitution.jl")
include("site.jl")
include("sequence_simulator.jl")
include("SequencePropagator.jl")
include("export.jl")


export JC, F81, K2P, HKY, GTR, SubstitutionModel
export SiteModel
export rand_seq, SequencePropagator
export Sequence
export write_fasta, write_nexus, write_phylip, write_alignment

end # module SeqSim
