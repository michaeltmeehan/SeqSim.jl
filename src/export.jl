using FilePathsBase  # or just use `splitext` manually

function format_label(seq::Sequence)
    label = isnothing(seq.taxon) ? "unknown" : string(seq.taxon)
    if !isnothing(seq.time)
        label *= "_$(round(seq.time; sigdigits=4))"
    end
    return label
end


function write_fasta(io::IO, alignment::Vector{Sequence})
    wrap_length = 60
    for seq in alignment
        println(io, ">", format_label(seq))
        dna = join(nucleotides[nt] for nt in seq.value)
        for i in 1:wrap_length:length(dna)
            println(io, dna[i:min(i+wrap_length-1, end)])
        end
    end
end


function write_nexus(io::IO, alignment::Vector{Sequence})
    n_taxa = length(alignment)
    seq_length = length(alignment[1].value)

    println(io, "#NEXUS\n")
    println(io, "BEGIN DATA;")
    println(io, "    DIMENSIONS NTAX=$n_taxa NCHAR=$seq_length;")
    println(io, "    FORMAT DATATYPE=DNA MISSING=? GAP=-;")
    println(io, "    MATRIX")

    for seq in alignment
        label = format_label(seq)
        dna = join(nucleotides[nt] for nt in seq.value)
        println(io, "    $label    $dna")
    end

    println(io, "    ;")
    println(io, "END;")
end


function write_phylip(io::IO, alignment::Vector{Sequence})
    n_taxa = length(alignment)
    seq_length = length(alignment[1].value)

    println(io, "$n_taxa $seq_length")

    for seq in alignment
        label = format_label(seq)
        short_label = length(label) > 10 ? first(label, 10) : rpad(label, 10)
        dna = join(nucleotides[nt] for nt in seq.value)
        println(io, "$short_label$dna")
    end
end

function infer_format_from_extension(filename::String)
    ext = lowercase(splitext(filename)[2])
    if ext in [".fasta", ".fa"]
        return :fasta
    elseif ext in [".nex", ".nexus"]
        return :nexus
    elseif ext in [".phy", ".phylip"]
        return :phylip
    else
        error("Cannot infer format from file extension: $ext")
    end
end


function write_alignment(filename::String, alignment::Vector{Sequence}; format::Union{Symbol, Nothing}=nothing)
    format = isnothing(format) ? infer_format_from_extension(filename) : format

    open(filename, "w") do io
        if format == :fasta
            write_fasta(io, alignment)
        elseif format == :nexus
            write_nexus(io, alignment)
        elseif format == :phylip
            write_phylip(io, alignment)
        else
            error("Unsupported format: $format")
        end
    end
end