@testset "sequence and alignment representation" begin
    seq = Sequence("ACGT"; taxon="sample-1", time=1.5)
    @test seq.value == "ACGT"
    @test seq.taxon == "sample-1"
    @test seq.time == 1.5
    @test Sequence(UInt8[1, 2, 3, 4]).value == "ACGT"
    @test SeqSim.encode("TGCA") == UInt8[4, 3, 2, 1]
    @test SeqSim.decode(UInt8[1, 2, 3, 4]) == "ACGT"

    @test_throws ArgumentError Sequence("")
    @test_throws ArgumentError Sequence("ACGN")
    @test_throws ArgumentError Sequence("ACGT"; time=Inf)
    @test_throws ArgumentError Sequence(UInt8[1, 5])
    @test_throws ArgumentError SeqSim.encode("AC-GT")
    @test_throws ArgumentError SeqSim.decode(UInt8[0])

    valid_alignment = [Sequence("ACGT"; taxon="a"), Sequence("AGGT"; taxon="b")]
    @test SeqSim.get_snps(valid_alignment) == [2]

    @test_throws ArgumentError SeqSim.get_snps(Sequence[])
    @test_throws ArgumentError SeqSim.get_snps([Sequence("ACGT"), Sequence("ACG")])
end

@testset "write-only alignment IO integrity" begin
    alignment = [
        Sequence("ACGTACGT"; taxon="root", time=0.0),
        Sequence("ACGTTCGT"; taxon="child", time=0.5),
    ]

    fasta = sprint(io -> write_fasta(io, alignment))
    @test occursin(">root_0.0", fasta)
    @test occursin(">child_0.5", fasta)
    @test occursin("ACGTACGT", fasta)
    @test occursin("ACGTTCGT", fasta)

    nexus = sprint(io -> write_nexus(io, alignment))
    @test occursin("#NEXUS", nexus)
    @test occursin("NTAX=2", nexus)
    @test occursin("NCHAR=8", nexus)
    @test occursin("root_0.0    ACGTACGT", nexus)

    phylip = sprint(io -> write_phylip(io, alignment))
    @test startswith(phylip, "2 8\n")
    @test occursin("root_0.0  ACGTACGT", phylip)

    @test_throws ArgumentError write_fasta(IOBuffer(), Sequence[])
    @test_throws ArgumentError write_nexus(IOBuffer(), [Sequence("ACGT"; taxon="a"), Sequence("AC"; taxon="b")])
    @test_throws ArgumentError write_phylip(IOBuffer(), [Sequence("ACGT"), Sequence("ACGT")])
    @test_throws ArgumentError write_phylip(IOBuffer(), [Sequence("ACGT"; taxon="abcdefghijX"), Sequence("ACGT"; taxon="abcdefghijY")])
end
