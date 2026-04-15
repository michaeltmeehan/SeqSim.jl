"""
The supported DNA alphabet, in the same order as the internal state encoding.

`SeqSim.jl` currently supports unambiguous DNA sequences only. The stable
external representation is an `AbstractString` containing only `A`, `C`, `G`,
and `T`; the internal simulator representation uses `UInt8(1):UInt8(4)`.
Ambiguity codes, gaps, RNA, amino acids, and missing data are intentionally
unsupported in the current standalone core.
"""
const nucleotides = ['A', 'C', 'G', 'T']
const nucleotide_map = Dict('A' => UInt8(1), 'C' => UInt8(2), 'G' => UInt8(3), 'T' => UInt8(4))

function decode(i::UInt8)
    i in UInt8(1):UInt8(4) || throw(ArgumentError("Encoded nucleotide state must be in 1:4; got $i."))
    return nucleotides[i]
end


function decode(v::Vector{UInt8})
    validate_encoded_sequence(v)
    return join(nucleotides[i] for i in v)
end


function encode(nucleotide::Char)
    haskey(nucleotide_map, nucleotide) || throw(ArgumentError("Unsupported DNA symbol '$nucleotide'. Supported symbols are A, C, G, and T."))
    return nucleotide_map[nucleotide]
end


function encode(sequence::AbstractString)
    validate_dna_sequence(sequence)
    return [encode(nucleotide) for nucleotide in sequence]
end


function validate_dna_sequence(sequence::AbstractString)
    isempty(sequence) && throw(ArgumentError("Sequence value must not be empty."))
    for nucleotide in sequence
        haskey(nucleotide_map, nucleotide) || throw(ArgumentError("Unsupported DNA symbol '$nucleotide'. Supported symbols are A, C, G, and T."))
    end
    return sequence
end


function validate_encoded_sequence(sequence::Vector{UInt8})
    isempty(sequence) && throw(ArgumentError("Encoded sequence must not be empty."))
    for state in sequence
        state in UInt8(1):UInt8(4) || throw(ArgumentError("Encoded nucleotide state must be in 1:4; got $state."))
    end
    return sequence
end


"""
A dictionary mapping nucleotide characters to their corresponding color representations.

# Mappings
- `'A'`: Green
- `'C'`: Blue
- `'G'`: Magenta (commonly used for purple)
- `'T'`: Red

Each nucleotide is associated with a `Crayon` object that specifies its foreground color for terminal output.
"""
nucleotide_colors = Dict(
    'A' => Crayon(foreground = :green),
    'C' => Crayon(foreground = :blue),
    'G' => Crayon(foreground = :magenta),  # Magenta is often used for purple
    'T' => Crayon(foreground = :red)
)


# Styles
title_style = Crayon(bold=true)  # Just bold for titles
value_style = Crayon(foreground=:cyan)
snp_style = Crayon(foreground=:yellow)
note_style = Crayon(foreground=:white)  # dimmer note for "more taxa not shown"


struct Sequence
    taxon::Union{Nothing, AbstractString, Int}
    value::AbstractString
    time::Union{Nothing, Float64}

    function Sequence(taxon::Union{Nothing, AbstractString, Int}, value::AbstractString, time::Union{Nothing, Float64})
        validate_dna_sequence(value)
        time !== nothing && isfinite(time) || time === nothing || throw(ArgumentError("Sequence time must be finite when provided."))
        return new(taxon, value, time)
    end
end


function Sequence(seq::AbstractString; taxon=nothing, time=nothing)
    return Sequence(taxon, seq, time)
end


function Sequence(seq::Vector{UInt8}; taxon=nothing, time=nothing)
    value = decode(seq)
    return Sequence(taxon, value, time)
end


const Alignment = Vector{Sequence}


function validate_alignment(alignment::Alignment)
    isempty(alignment) && throw(ArgumentError("Alignment must contain at least one sequence."))
    seq_length = length(first(alignment).value)
    for (index, seq) in pairs(alignment)
        length(seq.value) == seq_length || throw(ArgumentError("Alignment sequence $index has length $(length(seq.value)); expected $seq_length."))
        validate_dna_sequence(seq.value)
    end
    return alignment
end


function color_sequence(io::IO, sequence::AbstractString)
    for nucleotide in sequence
        crayon = get(nucleotide_colors, nucleotide, nothing)
        if crayon !== nothing
            print(io, crayon(string(nucleotide)))
        else
            print(io, nucleotide)
        end
    end
end


function showcompact(io, seq::AbstractString)
    if isempty(seq)
        println(io, "< EMPTY SEQUENCE >")
    else
        width = displaysize()[2]
        if length(seq) > width
            half = div(width, 2)
            color_sequence(io, seq[1:half-1])
            print(io, "...")
            color_sequence(io, seq[end-half+2:end])
        else
            color_sequence(io, seq)
        end
    end
end


function Base.show(io::IO, seq::Sequence)
    println(io, value_style("$(length(seq.value)) bp "), "nt sequence")
    showcompact(io, seq.value)    
    seq.taxon !== nothing && print(io, "\ntaxon=$(seq.taxon), ")
    seq.time !== nothing && print(io, "\ntime=$(round(seq.time, sigdigits=3))")
end


function Base.show(io::IO, ::MIME"text/plain", aln::Vector{Sequence})
    if isempty(aln)
        println(io, title_style("Alignment (empty)"))
        return
    end

    num_taxa = length(aln)
    seq_length = isempty(aln) ? 0 : length(aln[1].value)

    println(io, title_style("Alignment with "), value_style("$(num_taxa) taxa"))
    println(io, title_style("Sequence length : "), value_style("$seq_length"))
    
    snps = get_snps(aln)
    num_snps = length(snps)
    println(io, title_style("Number of SNPs  : "), value_style("$num_snps"))
    

    if num_snps > 0 
        println(io, title_style("SNP sites       : "), snp_style(join(first(snps, min(5, num_snps)), ", ")), 
        num_snps > 5 ? snp_style(", ...") : "")
    end
    println(io)

    # Thresholds for truncation
    max_taxa_to_display = 5

    truncate_taxa = num_taxa > max_taxa_to_display

    num_taxa_display = truncate_taxa ? max_taxa_to_display : num_taxa

    taxa = [seq.taxon === nothing ? "?" : string(seq.taxon) for seq in aln[1:num_taxa_display]]
    pad = maximum(length.(taxa)) + 1

    for seq in aln[1:num_taxa_display]
        label = isnothing(seq.taxon) ? "unknown" : string(seq.taxon)
        print(io, rpad(label, pad), ": ")
        showcompact(io, seq.value)
        println(io)
    end

    if truncate_taxa
        println(io, note_style("... (" * string(num_taxa - num_taxa_display) * " more taxa not shown)"))
    end
end
