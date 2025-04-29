"""
A constant array representing the four standard nucleotides in DNA: Adenine ('A'), Cytosine ('C'), Guanine ('G'), and Thymine ('T').
"""
const nucleotides = ['A', 'C', 'G', 'T']


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


"""
    Sequence

A mutable struct representing a biological sequence or similar data structure.

# Fields
- `taxon::Union{Nothing, String, Int}`: An identifier for the sequence, which can be `Nothing`, a `String`, or an `Int`.
- `value::Vector{UInt8}`: The sequence data, stored as a vector of `UInt8` values.
- `time::Union{Nothing, Float64}`: An optional timestamp or time-related value associated with the sequence, which can be `Nothing` or a `Float64`.
"""
mutable struct Sequence
    taxon::Union{Nothing, String, Int}
    value::Vector{UInt8}
    time::Union{Nothing, Float64}
end


"""
    Sequence(seq::Vector{UInt8}; taxon=nothing, time=nothing, metadata=Dict{Symbol,Any}())

Constructs a `Sequence` object.

# Arguments
- `seq::Vector{UInt8}`: The sequence data represented as a vector of `UInt8`.
- `taxon`: An optional identifier for the sequence. Defaults to `nothing`.
- `time`: An optional timestamp or time-related information for the sequence. Defaults to `nothing`.

# Returns
A `Sequence` object initialized with the provided sequence data and optional parameters.
"""
function Sequence(seq::Vector{UInt8}; taxon=nothing, time=nothing)
    return Sequence(taxon, seq, time)
end


"""
    color_sequence(io::IO, sequence::AbstractString)

Prints a color-coded representation of a nucleotide sequence to the given IO stream.

# Arguments
- `io::IO`: The IO stream where the colored sequence will be printed.
- `sequence::AbstractString`: The nucleotide sequence to be color-coded and printed.

# Details
Each nucleotide in the sequence is assigned a color based on the `nucleotide_colors` dictionary. If a nucleotide has a corresponding color (represented as a `Crayon` object), it is printed in that color. If no color is defined for a nucleotide, it is printed without any color.

"""
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


"""
    Base.show(io::IO, seq::Sequence)

Custom `show` method for the `Sequence` type. This method defines how a `Sequence` object
is displayed when printed to an output stream.

# Arguments
- `io::IO`: The output stream where the `Sequence` object will be printed.
- `seq::Sequence`: The `Sequence` object to be displayed.

# Behavior
- Displays the `Sequence` object in the format `Sequence(...)`.
- Includes the `taxon` field if it is not `nothing`.
- Displays the `value` field as a string of nucleotides, with optional coloring applied.
- Includes the `time` field if it is not `nothing`.

This method ensures a clear and informative representation of the `Sequence` object for debugging
or logging purposes.
"""
function Base.show(io::IO, seq::Sequence)
    print(io, "Sequence(")
    seq.taxon !== nothing && print(io, "taxon=$(seq.taxon), ")
    # print(io, "value=\"", join(nucleotides[nucl] for nucl in seq.value), "\"")
    print(io, "value=\"")
    color_sequence(io, join(nucleotides[nucl] for nucl in seq.value))
    print(io, "\"")
    seq.time !== nothing && print(io, ", time=$(round(seq.time, sigdigits=3))")
    print(io, ")")
end


function get_snps(aln::Vector{Sequence})
    snps = Vector{Int}()

    for site in eachindex(aln[1].value)
        nucleotide = aln[1].value[site]

        for seq in aln[2:end]
            if seq.value[site] != nucleotide
                push!(snps, site)
                break
            end
        end
    end
    return snps
end


"""
    Base.show(io::IO, ::MIME"text/plain", aln::Vector{Sequence})

Custom `show` method for displaying a vector of `Sequence` objects (an alignment) in a human-readable format.

# Arguments
- `io::IO`: The output stream where the alignment will be printed.
- `::MIME"text/plain"`: Specifies that the output is intended for plain text display.
- `aln::Vector{Sequence}`: A vector of `Sequence` objects representing the alignment.

# Behavior
- If the alignment is empty, it prints "Alignment (empty)".
- For non-empty alignments:
  - Each sequence is decoded into its nucleotide representation.
  - Sequence IDs are displayed alongside their corresponding nucleotide sequences.
  - IDs are padded to align the output for better readability.
  - Nucleotides are optionally colorized using a predefined mapping (`nucleotide_colors`).

# Notes
- Sequence IDs default to "?" if they are `nothing`.
- The method uses `rpad` to align sequence IDs and `crayon` for optional nucleotide colorization.
"""
# Styles
title_style = Crayon(bold=true)  # Just bold for titles
value_style = Crayon(foreground=:cyan)
snp_style = Crayon(foreground=:yellow)
note_style = Crayon(foreground=:white)  # dimmer note for "more taxa not shown"

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
    max_seq_length_to_display = 60

    truncate_taxa = num_taxa > max_taxa_to_display
    truncate_seq = seq_length > max_seq_length_to_display

    num_taxa_display = truncate_taxa ? max_taxa_to_display : num_taxa
    seq_length_display = truncate_seq ? max_seq_length_to_display : seq_length

    decoded = [join(nucleotides[nucl] for nucl in seq.value[1:seq_length_display]) for seq in aln[1:num_taxa_display]]
    taxa = [seq.taxon === nothing ? "?" : string(seq.taxon) for seq in aln[1:num_taxa_display]]
    pad = maximum(length.(taxa)) + 1

    for (taxon, seq_str) in zip(taxa, decoded)
        print(io, rpad(taxon, pad), ": ")
        for nt in seq_str
            crayon = get(nucleotide_colors, nt, identity)
            print(io, crayon(string(nt)))
        end
        if truncate_seq
            print(io, " ...")
        end
        println(io)
    end

    if truncate_taxa
        println(io, note_style("... (" * string(num_taxa - num_taxa_display) * " more taxa not shown)"))
    end
end
