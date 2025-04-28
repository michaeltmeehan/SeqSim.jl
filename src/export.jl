using FilePathsBase  # or just use `splitext` manually

"""
    write_fasta(filename::String, alignment::Vector{Sequence})

Writes an alignment to a FASTA file.

If a sequence has a non-`nothing` `time` field, it appends it to the taxon label separated by an underscore.
"""
function write_fasta(filename::String, alignment::Vector{Sequence})
    open(filename, "w") do io
        for seq in alignment
            # Construct taxon label
            label = isnothing(seq.taxon) ? "unknown" : string(seq.taxon)
            if !isnothing(seq.time)
                label *= "_" * string(seq.time)
            end
            println(io, ">", label)

            # Build full sequence string
            seq_str = join(nucleotides[nucl] for nucl in seq.value)

            # Write the sequence with wrapping at 60 characters
            wrap_length = 60
            for i in 1:wrap_length:length(seq_str)
                println(io, seq_str[i:min(i+wrap_length-1, end)])
            end
        end
    end
end


function write_nexus(filename::String, alignment::Vector{Sequence})
    n_taxa = length(alignment)
    seq_length = length(alignment[1].value)

    open(filename, "w") do io
        # NEXUS header
        println(io, "#NEXUS")
        println(io)
        println(io, "BEGIN DATA;")
        println(io, "    DIMENSIONS NTAX=$n_taxa NCHAR=$seq_length;")
        println(io, "    FORMAT DATATYPE=DNA MISSING=? GAP=-;")
        println(io, "    MATRIX")

        # Sequences
        for seq in alignment
            # Build label
            label = isnothing(seq.taxon) ? "unknown" : string(seq.taxon)
            if !isnothing(seq.time)
                label *= "_" * string(seq.time)
            end

            # Build full sequence string
            seq_str = join(nucleotides[nucl] for nucl in seq.value)

            # Write label and sequence
            println(io, label, " ", seq_str)
        end

        # End block
        println(io, "    ;")
        println(io, "END;")
    end
end


function write_phylip(filename::String, alignment::Vector{Sequence})
    n_taxa = length(alignment)
    seq_length = length(alignment[1].value)

    open(filename, "w") do io
        # PHYLIP header
        println(io, "$n_taxa $seq_length")

        # Sequences
        for seq in alignment
            # Build label
            label = isnothing(seq.taxon) ? "unknown" : string(seq.taxon)
            if !isnothing(seq.time)
                label *= "_" * string(seq.time)
            end

            # PHYLIP format requires taxon labels <=10 characters
            short_label = length(label) > 10 ? first(label, 10) : rpad(label, 10)

            # Build full sequence string
            seq_str = join(nucleotides[nucl] for nucl in seq.value)

            println(io, short_label * seq_str)
        end
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
    if format === nothing
        format = infer_format_from_extension(filename)
    end

    if format == :fasta
        write_fasta(filename, alignment)
    elseif format == :nexus
        write_nexus(filename, alignment)
    elseif format == :phylip
        write_phylip(filename, alignment)
    else
        error("Unsupported format: $format. Supported formats are :fasta, :nexus, :phylip.")
    end
end
