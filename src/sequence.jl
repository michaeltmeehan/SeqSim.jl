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
- `id::Union{Nothing, String, Int}`: An identifier for the sequence, which can be `Nothing`, a `String`, or an `Int`.
- `sequence::Vector{UInt8}`: The sequence data, stored as a vector of `UInt8` values.
- `time::Union{Nothing, Float64}`: An optional timestamp or time-related value associated with the sequence, which can be `Nothing` or a `Float64`.
- `metadata::Dict{Symbol, Any}`: A dictionary for storing additional metadata about the sequence, where keys are `Symbol`s and values can be of any type.
"""
mutable struct Sequence
    id::Union{Nothing, String, Int}
    sequence::Vector{UInt8}
    time::Union{Nothing, Float64}
    metadata::Dict{Symbol, Any}
end


"""
    Sequence(seq::Vector{UInt8}; id=nothing, time=nothing, metadata=Dict{Symbol,Any}())

Constructs a `Sequence` object.

# Arguments
- `seq::Vector{UInt8}`: The sequence data represented as a vector of `UInt8`.
- `id`: An optional identifier for the sequence. Defaults to `nothing`.
- `time`: An optional timestamp or time-related information for the sequence. Defaults to `nothing`.
- `metadata::Dict{Symbol,Any}`: An optional dictionary containing metadata associated with the sequence. Defaults to an empty dictionary.

# Returns
A `Sequence` object initialized with the provided sequence data and optional parameters.
"""
function Sequence(seq::Vector{UInt8}; id=nothing, time=nothing, metadata=Dict{Symbol,Any}())
    return Sequence(id, seq, time, metadata)
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
- Includes the `id` field if it is not `nothing`.
- Displays the `sequence` field as a string of nucleotides, with optional coloring applied.
- Includes the `time` field if it is not `nothing`.
- Includes the `metadata` field if it is not empty.

This method ensures a clear and informative representation of the `Sequence` object for debugging
or logging purposes.
"""
function Base.show(io::IO, seq::Sequence)
    print(io, "Sequence(")
    seq.id !== nothing && print(io, "id=$(seq.id), ")
    # print(io, "sequence=\"", join(nucleotides[nucl] for nucl in seq.sequence), "\"")
    print(io, "sequence=\"")
    color_sequence(io, join(nucleotides[nucl] for nucl in seq.sequence))
    print(io, "\"")
    seq.time !== nothing && print(io, ", time=$(seq.time)")
    !isempty(seq.metadata) && print(io, ", metadata=$(seq.metadata)")
    print(io, ")")
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
function Base.show(io::IO, ::MIME"text/plain", aln::Vector{Sequence})
    if isempty(aln)
        println(io, "Alignment (empty)")
        return
    end

    # Decode sequences to nucleotides
    decoded = [join(nucleotides[nucl] for nucl in s.sequence) for s in aln]
    ids = [s.id === nothing ? "?" : string(s.id) for s in aln]
    pad = maximum(length.(ids))

    for (id, seq_str) in zip(ids, decoded)
        print(io, rpad(id, pad), ": ")
        for nt in seq_str
            crayon = get(nucleotide_colors, nt, identity)  # fallback: identity(x) = x
            print(io, crayon(string(nt)))
        end
        println(io)
    end
end



