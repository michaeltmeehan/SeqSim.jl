"""
A constant array representing the four standard nucleotides in DNA: Adenine ('A'), Cytosine ('C'), Guanine ('G'), and Thymine ('T').
"""
const nucleotides = ['A', 'C', 'G', 'T']
const nucleotide_map = Dict(0x01 => 'A', 0x02 => 'C', 0x03 =>  'G', 0x04 => 'T')

function decode(i::UInt8)
    return nucleotides[i]
end


function decode(v::Vector{UInt8})
    return join(nucleotides[i] for i in v)
end


function encode(nucleotide::Char)
    return nucleotide_map[nucleotide]
end


function encode(sequence::AbstractString)
    return [encode(nucleotide) for nucleotide in sequence]
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
end


function Sequence(seq::AbstractString; taxon=nothing, time=nothing)
    return Sequence(taxon, seq, time)
end


function Sequence(seq::Vector{UInt8}; taxon=nothing, time=nothing)
    value = decode(seq)
    return Sequence(taxon, value, time)
end


const Alignment = Vector{Sequence}


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


function get_snps(aln::Alignment)
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
