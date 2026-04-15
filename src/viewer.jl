struct AlignmentViewer
    alignment::Alignment
    inspection::AlignmentInspection
    reference::Any
    max_sequences::Int
    max_sites::Int
    variable_only::Bool
    show_site_summary::Bool
end


function validate_preview_limit(name::AbstractString, value::Integer)
    value > 0 || throw(ArgumentError("$name must be positive; got $value."))
    return Int(value)
end


function AlignmentViewer(
    alignment::Alignment;
    sites::SiteSelector=:,
    reference=nothing,
    max_sequences::Integer=10,
    max_sites::Integer=80,
    variable_only::Bool=false,
    show_site_summary::Bool=true,
)
    validate_alignment(alignment)
    inspection = inspect_sites(alignment; sites=sites, reference=reference)
    return AlignmentViewer(
        alignment,
        inspection,
        reference,
        validate_preview_limit("max_sequences", max_sequences),
        validate_preview_limit("max_sites", max_sites),
        variable_only,
        show_site_summary,
    )
end


alignment_viewer(alignment::Alignment; kwargs...) = AlignmentViewer(alignment; kwargs...)


function window_viewer(alignment::Alignment, first_site::Integer, last_site::Integer; kwargs...)
    return AlignmentViewer(alignment; sites=site_window(alignment, first_site, last_site), kwargs...)
end


function centered_window_viewer(alignment::Alignment, focal_site::Integer; radius::Integer=5, kwargs...)
    return AlignmentViewer(alignment; sites=centered_site_window(alignment, focal_site; radius=radius), kwargs...)
end


function variable_site_viewer(alignment::Alignment; sites::SiteSelector=:, kwargs...)
    selected_sites = normalize_site_selector(sites, length(first(validate_alignment(alignment)).value); allow_empty=true)
    variable = Set(variable_sites(alignment))
    return AlignmentViewer(alignment; sites=[site for site in selected_sites if site in variable], kwargs...)
end


function visible_rows(viewer::AlignmentViewer)
    if viewer.variable_only
        return [row for row in viewer.inspection.rows if row.is_variable]
    end
    return viewer.inspection.rows
end


visible_sites(viewer::AlignmentViewer) = [row.site for row in visible_rows(viewer)]


function displayed_sites(viewer::AlignmentViewer)
    sites = visible_sites(viewer)
    return sites[1:min(length(sites), viewer.max_sites)]
end


function reference_string(viewer::AlignmentViewer, sites::Vector{Int})
    viewer.reference === nothing && return nothing
    return filtered_value(reference_value(viewer.alignment, viewer.reference), sites)
end


function viewer_preview(viewer::AlignmentViewer)
    return sprint(io -> show_viewer_preview(io, viewer))
end


function show_viewer_preview(io::IO, viewer::AlignmentViewer)
    rows = visible_rows(viewer)
    sites = [row.site for row in rows]
    shown_sites = sites[1:min(length(sites), viewer.max_sites)]
    shown_rows = rows[1:min(length(rows), viewer.max_sites)]
    displayed_sequence_count = min(length(viewer.alignment), viewer.max_sequences)

    println(io, "AlignmentViewer")
    println(io, "  Sequences      : ", length(viewer.alignment), " (showing ", displayed_sequence_count, ")")
    println(io, "  Selected sites : ", length(sites), " (showing ", length(shown_sites), ")")
    viewer.variable_only && println(io, "  Site filter    : variable only")
    viewer.inspection.reference !== nothing && println(io, "  Reference      : ", viewer.inspection.reference)

    if isempty(shown_sites)
        println(io, "  Sites          : <none>")
        return
    end

    println(io, "Sites: ", join(shown_sites, " "))

    if viewer.show_site_summary
        maf = [round(row.minor_allele_frequency; sigdigits=3) for row in shown_rows]
        variable = [row.is_variable ? "*" : "." for row in shown_rows]
        println(io, "Variable: ", join(variable, " "))
        println(io, "MAF: ", join(maf, " "))
        if viewer.reference !== nothing
            nonref = [round(row.nonreference_frequency; sigdigits=3) for row in shown_rows]
            println(io, "Nonref: ", join(nonref, " "))
        end
    end

    ref = reference_string(viewer, shown_sites)
    ref !== nothing && println(io, rpad("reference", 12), ref)

    projections = selected_site_strings(viewer.alignment; sites=shown_sites)
    for projection in projections[1:displayed_sequence_count]
        label = projection.taxon === nothing ? "unknown" : string(projection.taxon)
        println(io, rpad(label, 12), projection.value)
        if ref !== nothing
            markers = String([projection.value[i] == ref[i] ? UInt8('.') : UInt8('^') for i in eachindex(projection.value)])
            println(io, rpad("", 12), markers)
        end
    end

    length(viewer.alignment) > displayed_sequence_count &&
        println(io, "... (", length(viewer.alignment) - displayed_sequence_count, " more sequences)")
    length(sites) > length(shown_sites) &&
        println(io, "... (", length(sites) - length(shown_sites), " more selected sites)")
end


function Base.show(io::IO, viewer::AlignmentViewer)
    print(io, "AlignmentViewer(",
        length(viewer.alignment), " sequences, ",
        length(visible_sites(viewer)), " selected sites")
    viewer.variable_only && print(io, ", variable_only=true")
    viewer.inspection.reference !== nothing && print(io, ", reference=", viewer.inspection.reference)
    print(io, ")")
end


Base.show(io::IO, ::MIME"text/plain", viewer::AlignmentViewer) = show_viewer_preview(io, viewer)
