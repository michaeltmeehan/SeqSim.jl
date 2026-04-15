const SiteSelector = Union{Colon, AbstractVector{<:Integer}, AbstractRange{<:Integer}, AbstractVector{Bool}}

"""
Computation-facing DNA state order used by SeqSim processing utilities.

The supported DNA alphabet is encoded as `A => 1`, `C => 2`, `G => 3`, and
`T => 4`. Functions such as `site_states` and `site_state_counts` use this
state order so downstream analysis can avoid display-oriented character work.
"""
const DNA_STATE_ORDER = (:A, :C, :G, :T)

state_char(state::Integer) = decode(UInt8(state))
state_index(nucleotide::Char) = Int(encode(nucleotide))


function validate_site_indices(sites::AbstractVector{<:Integer}, sequence_length::Integer; allow_empty::Bool=false)
    isempty(sites) && !allow_empty && throw(ArgumentError("At least one site index must be provided."))

    seen = Set{Int}()
    validated = Vector{Int}(undef, length(sites))
    for (i, site) in pairs(sites)
        1 <= site <= sequence_length || throw(ArgumentError("Site index $site is outside the valid range 1:$sequence_length."))
        site_int = Int(site)
        site_int in seen && throw(ArgumentError("Site indices must be unique; duplicate site $site_int."))
        push!(seen, site_int)
        validated[i] = site_int
    end
    return validated
end


validate_site_indices(sites::AbstractRange{<:Integer}, sequence_length::Integer; allow_empty::Bool=false) =
    validate_site_indices(collect(sites), sequence_length; allow_empty=allow_empty)


function normalize_site_selector(sites::Colon, sequence_length::Integer; allow_empty::Bool=false)
    return collect(1:Int(sequence_length))
end


function normalize_site_selector(sites::AbstractVector{Bool}, sequence_length::Integer; allow_empty::Bool=false)
    length(sites) == sequence_length || throw(ArgumentError("Boolean site mask length must be $sequence_length; got $(length(sites))."))
    selected = findall(sites)
    isempty(selected) && !allow_empty && throw(ArgumentError("At least one site index must be provided."))
    return selected
end


normalize_site_selector(sites::AbstractVector{<:Integer}, sequence_length::Integer; allow_empty::Bool=false) =
    validate_site_indices(sites, sequence_length; allow_empty=allow_empty)


normalize_site_selector(sites::AbstractRange{<:Integer}, sequence_length::Integer; allow_empty::Bool=false) =
    validate_site_indices(sites, sequence_length; allow_empty=allow_empty)


function normalize_site_selector(sites, sequence_length::Integer; allow_empty::Bool=false)
    throw(ArgumentError("Unsupported site selector type $(typeof(sites)). Use :, an integer vector or range, or a boolean mask."))
end


function filter_sequences(predicate::Function, alignment::Alignment)
    validate_alignment(alignment)
    return [seq for seq in alignment if predicate(seq)]
end


function filtered_value(sequence::AbstractString, sites::AbstractVector{Int})
    bytes = Vector{UInt8}(undef, length(sites))
    for (i, site) in pairs(sites)
        bytes[i] = codeunit(sequence, site)
    end
    return String(bytes)
end


function filter_sites(alignment::Alignment, sites::SiteSelector)
    validate_alignment(alignment)
    selected_sites = normalize_site_selector(sites, length(first(alignment).value))
    return [
        Sequence(seq.taxon, filtered_value(seq.value, selected_sites), seq.time)
        for seq in alignment
    ]
end


"""
    site_states(alignment, site)

Return the encoded DNA states observed at `site` as `Vector{UInt8}` using the
internal `1:4` state order for `A`, `C`, `G`, and `T`.
"""
function site_states(alignment::Alignment, site::Integer)
    validate_alignment(alignment)
    selected_site = only(validate_site_indices([site], length(first(alignment).value)))
    states = Vector{UInt8}(undef, length(alignment))
    for (i, seq) in pairs(alignment)
        states[i] = encode(seq.value[selected_site])
    end
    return states
end


function site_state_counts(alignment::Alignment; sites::SiteSelector=:)
    validate_alignment(alignment)
    selected_sites = normalize_site_selector(sites, length(first(alignment).value); allow_empty=true)
    counts = zeros(Int, 4, length(selected_sites))
    for seq in alignment
        value = seq.value
        for (j, site) in pairs(selected_sites)
            counts[state_index(value[site]), j] += 1
        end
    end
    return counts
end


function observed_state_counts(alignment::Alignment; sites::SiteSelector=:)
    counts = site_state_counts(alignment; sites=sites)
    return [counts[:, site] for site in axes(counts, 2)]
end


function variable_sites(alignment::Alignment)
    counts = site_state_counts(alignment)
    return [site for site in axes(counts, 2) if count(!iszero, view(counts, :, site)) > 1]
end


function invariant_sites(alignment::Alignment)
    counts = site_state_counts(alignment)
    return [site for site in axes(counts, 2) if count(!iszero, view(counts, :, site)) == 1]
end


"""
    get_snps(alignment)

Compatibility alias for `variable_sites(alignment)`.

With the current unambiguous DNA alphabet, SNP sites and variable sites are the
same. Future gap or ambiguity support should extend the variable-site
classification first, then decide which variable sites qualify as SNPs.
"""
get_snps(alignment::Alignment) = variable_sites(alignment)


function snp_alignment(alignment::Alignment)
    sites = variable_sites(alignment)
    isempty(sites) && throw(ArgumentError("Alignment has no variable/SNP sites to extract."))
    return filter_sites(alignment, sites)
end


invariant_filtered_alignment(alignment::Alignment) = snp_alignment(alignment)


function reference_value(alignment::Alignment, reference::Integer)
    1 <= reference <= length(alignment) || throw(ArgumentError("Reference sequence index must be in 1:$(length(alignment)); got $reference."))
    return alignment[reference].value
end


function reference_value(alignment::Alignment, reference::Sequence)
    length(reference.value) == length(first(alignment).value) || throw(ArgumentError("Reference sequence length must be $(length(first(alignment).value)); got $(length(reference.value))."))
    return reference.value
end


function reference_value(alignment::Alignment, reference::AbstractString)
    validate_dna_sequence(reference)
    length(reference) == length(first(alignment).value) || throw(ArgumentError("Reference sequence length must be $(length(first(alignment).value)); got $(length(reference))."))
    return reference
end


function reference_value(alignment::Alignment, reference)
    throw(ArgumentError("Unsupported reference type $(typeof(reference)). Use a sequence index, Sequence, or DNA string."))
end


function site_snp_counts(alignment::Alignment; reference::Union{Integer, Sequence, AbstractString}=1, sites::SiteSelector=:)
    validate_alignment(alignment)
    selected_sites = normalize_site_selector(sites, length(first(alignment).value); allow_empty=true)
    ref = reference_value(alignment, reference)
    counts = zeros(Int, length(selected_sites))
    for seq in alignment
        value = seq.value
        for (j, site) in pairs(selected_sites)
            counts[j] += value[site] != ref[site]
        end
    end
    return counts
end


function sequence_snp_counts(alignment::Alignment; reference::Union{Integer, Sequence, AbstractString}=1, sites::SiteSelector=:)
    validate_alignment(alignment)
    selected_sites = normalize_site_selector(sites, length(first(alignment).value); allow_empty=true)
    ref = reference_value(alignment, reference)
    counts = zeros(Int, length(alignment))
    for (i, seq) in pairs(alignment)
        value = seq.value
        for site in selected_sites
            counts[i] += value[site] != ref[site]
        end
    end
    return counts
end


function minor_allele_counts(alignment::Alignment; sites::SiteSelector=:)
    counts = site_state_counts(alignment; sites=sites)
    result = zeros(Int, size(counts, 2))
    for site in axes(counts, 2)
        observed = 0
        minor = typemax(Int)
        for count in view(counts, :, site)
            if count > 0
                observed += 1
                minor = min(minor, count)
            end
        end
        result[site] = observed <= 1 ? 0 : minor
    end
    return result
end


minor_allele_frequencies(alignment::Alignment; sites::SiteSelector=:) =
    minor_allele_counts(alignment; sites=sites) ./ length(validate_alignment(alignment))


function nonreference_frequencies(alignment::Alignment; reference::Union{Integer, Sequence, AbstractString}=1, sites::SiteSelector=:)
    return site_snp_counts(alignment; reference=reference, sites=sites) ./ length(validate_alignment(alignment))
end


function singleton_sites(alignment::Alignment)
    counts = site_state_counts(alignment)
    return [site for site in axes(counts, 2) if count(!iszero, view(counts, :, site)) > 1 && any(==(1), view(counts, :, site))]
end


function parsimony_informative_sites(alignment::Alignment)
    counts = site_state_counts(alignment)
    return [site for site in axes(counts, 2) if count(count -> count >= 2, view(counts, :, site)) >= 2]
end


function pairwise_differences(seq1::Sequence, seq2::Sequence; sites::SiteSelector=:)
    length(seq1.value) == length(seq2.value) || throw(ArgumentError("Sequences must have equal lengths; got $(length(seq1.value)) and $(length(seq2.value))."))
    selected_sites = normalize_site_selector(sites, length(seq1.value); allow_empty=true)
    differences = 0
    for site in selected_sites
        differences += seq1.value[site] != seq2.value[site]
    end
    return differences
end


function pairwise_difference_matrix(alignment::Alignment; sites::SiteSelector=:)
    validate_alignment(alignment)
    selected_sites = normalize_site_selector(sites, length(first(alignment).value); allow_empty=true)
    n = length(alignment)
    matrix = zeros(Int, n, n)
    for i in 1:n
        for j in i+1:n
            differences = 0
            value_i = alignment[i].value
            value_j = alignment[j].value
            for site in selected_sites
                differences += value_i[site] != value_j[site]
            end
            matrix[i, j] = differences
            matrix[j, i] = differences
        end
    end
    return matrix
end


function consensus_sequence(alignment::Alignment; taxon="consensus", time=nothing, tie::Symbol=:first)
    tie == :first || tie == :error || throw(ArgumentError("Unsupported consensus tie rule $tie. Use :first or :error."))
    counts = site_state_counts(alignment)
    sequence = Vector{UInt8}(undef, size(counts, 2))
    for site in axes(counts, 2)
        column = view(counts, :, site)
        max_count = maximum(column)
        winners = findall(==(max_count), column)
        if length(winners) > 1 && tie == :error
            throw(ArgumentError("Consensus tie at site $site."))
        end
        sequence[site] = UInt8(codeunit(string(state_char(first(winners))), 1))
    end
    return Sequence(taxon, String(sequence), time)
end


struct AlignmentSummary
    sequence_count::Int
    site_count::Int
    variable_site_count::Int
    invariant_site_count::Int
    singleton_site_count::Int
    parsimony_informative_site_count::Int
end


function AlignmentSummary(alignment::Alignment)
    validate_alignment(alignment)
    return AlignmentSummary(
        length(alignment),
        length(first(alignment).value),
        length(variable_sites(alignment)),
        length(invariant_sites(alignment)),
        length(singleton_sites(alignment)),
        length(parsimony_informative_sites(alignment)),
    )
end
