function site_window(alignment::Alignment, first_site::Integer, last_site::Integer)
    validate_alignment(alignment)
    sequence_length = length(first(alignment).value)
    1 <= first_site <= sequence_length || throw(ArgumentError("Window start must be in 1:$sequence_length; got $first_site."))
    1 <= last_site <= sequence_length || throw(ArgumentError("Window stop must be in 1:$sequence_length; got $last_site."))
    first_site <= last_site || throw(ArgumentError("Window start must be less than or equal to stop; got $first_site:$last_site."))
    return collect(Int(first_site):Int(last_site))
end


function centered_site_window(alignment::Alignment, focal_site::Integer; radius::Integer=5)
    validate_alignment(alignment)
    sequence_length = length(first(alignment).value)
    1 <= focal_site <= sequence_length || throw(ArgumentError("Focal site must be in 1:$sequence_length; got $focal_site."))
    radius >= 0 || throw(ArgumentError("Window radius must be non-negative; got $radius."))
    first_site = max(1, Int(focal_site) - Int(radius))
    last_site = min(sequence_length, Int(focal_site) + Int(radius))
    return collect(first_site:last_site)
end


struct SiteInspectionRow
    site::Int
    observed_states::Vector{UInt8}
    state_counts::Vector{Int}
    is_variable::Bool
    minor_allele_count::Int
    minor_allele_frequency::Float64
    is_singleton::Bool
    is_parsimony_informative::Bool
    nonreference_frequency::Union{Nothing, Float64}
end


struct AlignmentInspection
    sequence_count::Int
    site_count::Int
    selected_sites::Vector{Int}
    rows::Vector{SiteInspectionRow}
    reference::Union{Nothing, Int, String}
end


function reference_label(reference)
    reference === nothing && return nothing
    reference isa Integer && return Int(reference)
    reference isa Sequence && return isnothing(reference.taxon) ? "external" : string(reference.taxon)
    reference isa AbstractString && return "external"
    return string(typeof(reference))
end


"""
    inspect_sites(alignment; sites=:, reference=nothing)

Build a derived inspection object for selected original alignment sites.

This is an inspection/view-model helper, not an alternate alignment
representation. `selected_sites` and each `SiteInspectionRow.site` retain
original alignment coordinates, including caller-specified order.
"""
function inspect_sites(alignment::Alignment; sites::SiteSelector=:, reference=nothing)
    validate_alignment(alignment)
    sequence_count = length(alignment)
    site_count = length(first(alignment).value)
    selected_sites = normalize_site_selector(sites, site_count; allow_empty=true)

    counts = site_state_counts(alignment; sites=selected_sites)
    nonref_frequencies = reference === nothing ? nothing : nonreference_frequencies(alignment; reference=reference, sites=selected_sites)

    rows = Vector{SiteInspectionRow}(undef, length(selected_sites))
    for (j, site) in pairs(selected_sites)
        column = counts[:, j]
        observed = UInt8[state for state in UInt8(1):UInt8(4) if column[Int(state)] > 0]
        observed_count = count(!iszero, column)
        minor_count = 0
        if observed_count > 1
            minor_count = minimum(count for count in column if count > 0)
        end
        rows[j] = SiteInspectionRow(
            site,
            observed,
            collect(column),
            observed_count > 1,
            minor_count,
            minor_count / sequence_count,
            observed_count > 1 && any(==(1), column),
            count(count -> count >= 2, column) >= 2,
            nonref_frequencies === nothing ? nothing : nonref_frequencies[j],
        )
    end

    return AlignmentInspection(sequence_count, site_count, selected_sites, rows, reference_label(reference))
end


function inspect_window(alignment::Alignment, first_site::Integer, last_site::Integer; reference=nothing)
    return inspect_sites(alignment; sites=site_window(alignment, first_site, last_site), reference=reference)
end


function inspect_centered_window(alignment::Alignment, focal_site::Integer; radius::Integer=5, reference=nothing)
    return inspect_sites(alignment; sites=centered_site_window(alignment, focal_site; radius=radius), reference=reference)
end


"""
    selected_site_strings(alignment; sites=:)

Return derived per-sequence strings for selected original alignment sites.

This helper is intended for inspection and viewer rendering. It preserves
sequence metadata and selected-site order, but it is not a core alignment
representation and should not be used as a replacement for `Vector{Sequence}`.
"""
function selected_site_strings(alignment::Alignment; sites::SiteSelector=:)
    validate_alignment(alignment)
    selected_sites = normalize_site_selector(sites, length(first(alignment).value); allow_empty=true)
    return [
        (taxon = seq.taxon, value = filtered_value(seq.value, selected_sites))
        for seq in alignment
    ]
end


function Base.show(io::IO, inspection::AlignmentInspection)
    print(io, "AlignmentInspection(",
        inspection.sequence_count, " sequences, ",
        inspection.site_count, " sites, ",
        length(inspection.selected_sites), " selected sites")
    inspection.reference !== nothing && print(io, ", reference=", inspection.reference)
    print(io, ")")
end


function Base.show(io::IO, ::MIME"text/plain", inspection::AlignmentInspection)
    println(io, "AlignmentInspection")
    println(io, "  Sequences      : ", inspection.sequence_count)
    println(io, "  Source sites   : ", inspection.site_count)
    println(io, "  Selected sites : ", length(inspection.selected_sites))
    inspection.reference !== nothing && println(io, "  Reference      : ", inspection.reference)
    isempty(inspection.selected_sites) && return

    max_rows = min(5, length(inspection.rows))
    println(io, "  Site rows      :")
    for row in inspection.rows[1:max_rows]
        states = join(state_char.(row.observed_states), "")
        nonref = row.nonreference_frequency === nothing ? "" : ", nonref=$(round(row.nonreference_frequency; sigdigits=3))"
        println(io, "    site ", row.site,
            " states=", states,
            " variable=", row.is_variable,
            " maf=", round(row.minor_allele_frequency; sigdigits=3),
            nonref)
    end
    length(inspection.rows) > max_rows && println(io, "    ... (", length(inspection.rows) - max_rows, " more sites)")
end
