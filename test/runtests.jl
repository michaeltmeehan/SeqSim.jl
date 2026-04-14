using Random
using SeqSim
using Test

@testset "SeqSim smoke tests" begin
    rng = MersenneTwister(1)

    substitution_model = HKY([0.25, 0.25, 0.25, 0.25], 2.0)
    site_model = SiteModel(12, 0.1, 2, 1.0, 0.0, substitution_model)

    sequence = rand_seq(rng, site_model; taxon="root", time=0.0)
    @test sequence isa Sequence
    @test sequence.value isa AbstractString
    @test length(sequence.value) == site_model.sequence_length
    @test all(base in ['A', 'C', 'G', 'T'] for base in sequence.value)

    propagator = SequencePropagator(site_model)
    encoded = SeqSim.encode(sequence.value)
    propagated = propagator(rng, encoded, 0.5)
    child = Sequence(propagated; taxon="child", time=0.5)

    @test propagated isa Vector{UInt8}
    @test length(propagated) == site_model.sequence_length
    @test all(state in UInt8(1):UInt8(4) for state in propagated)
    @test length(child.value) == site_model.sequence_length

    alignment = [sequence, child]

    fasta = sprint(io -> write_fasta(io, alignment))
    @test occursin(">root_0.0", fasta)
    @test occursin(">child_0.5", fasta)
    @test occursin(sequence.value, fasta)
    @test occursin(child.value, fasta)

    nexus = sprint(io -> write_nexus(io, alignment))
    @test occursin("#NEXUS", nexus)
    @test occursin("NTAX=2", nexus)
    @test occursin("NCHAR=12", nexus)

    phylip = sprint(io -> write_phylip(io, alignment))
    @test startswith(phylip, "2 12")
end
