function format_label(seq::Sequence)
    label = isnothing(seq.taxon) ? "unknown" : string(seq.taxon)
    if !isnothing(seq.time)
        label *= "_$(round(seq.time; sigdigits=4))"
    end
    return label
end

dna_string(seq::Sequence) = seq.value

phylip_label(label::AbstractString) = length(label) > 10 ? first(label, 10) : rpad(label, 10)


function validate_export_alignment(alignment::Vector{Sequence})
    validate_alignment(alignment)
    labels = format_label.(alignment)
    length(unique(labels)) == length(labels) || throw(ArgumentError("Formatted alignment labels must be unique for export."))
    return alignment
end


function validate_phylip_alignment(alignment::Vector{Sequence})
    validate_export_alignment(alignment)
    labels = phylip_label.(format_label.(alignment))
    length(unique(labels)) == length(labels) || throw(ArgumentError("PHYLIP labels must be unique after 10-character truncation."))
    return alignment
end


function write_fasta(io::IO, alignment::Vector{Sequence})
    validate_export_alignment(alignment)
    wrap_length = 60
    for seq in alignment
        println(io, ">", format_label(seq))
        dna = dna_string(seq)
        for i in 1:wrap_length:length(dna)
            println(io, dna[i:min(i+wrap_length-1, end)])
        end
    end
end


function write_nexus(io::IO, alignment::Vector{Sequence})
    validate_export_alignment(alignment)
    n_taxa = length(alignment)
    seq_length = length(alignment[1].value)

    println(io, "#NEXUS\n")
    println(io, "BEGIN DATA;")
    println(io, "    DIMENSIONS NTAX=$n_taxa NCHAR=$seq_length;")
    println(io, "    FORMAT DATATYPE=DNA MISSING=? GAP=-;")
    println(io, "    MATRIX")

    for seq in alignment
        label = format_label(seq)
        dna = dna_string(seq)
        println(io, "    $label    $dna")
    end

    println(io, "    ;")
    println(io, "END;")
end


function write_phylip(io::IO, alignment::Vector{Sequence})
    validate_phylip_alignment(alignment)
    n_taxa = length(alignment)
    seq_length = length(alignment[1].value)

    println(io, "$n_taxa $seq_length")

    for seq in alignment
        label = format_label(seq)
        short_label = phylip_label(label)
        dna = dna_string(seq)
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
