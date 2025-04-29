struct SequencePropagator{N}
    site_model::SiteModel
    decomposition::RateMatrixDecomposition{N}
    transition_weights::Vector{Matrix{Float64}}
end


function SequencePropagator(site_model::SiteModel)
    Q = rate_matrix(site_model.substitution_model)
    N = size(Q, 1)
    decomposition = decompose(Q)
    transition_weights = [zeros(N, N) for _ in 1:site_model.gamma_category_count]
    return SequencePropagator{N}(site_model, decomposition, transition_weights)
end


function Base.show(io::IO, prop::SequencePropagator{N}) where {N}
    sm = prop.site_model
    println(io, "SequencePropagator{$N}")
    println(io, "  Sequence length     : ", sm.sequence_length)
    println(io, "  Mutation rate       : ", sm.mutation_rate)
    println(io, "  Substitution model  : ", typeof(sm.substitution_model))
    println(io, "  Gamma categories    : ", sm.gamma_category_count)
    println(io, "  Proportion invariant: ", sm.proportion_invariant)
end



function (prop::SequencePropagator{N})(rng::AbstractRNG, sequence::Vector{UInt8}, Δt::Float64) where {N}
    dec = prop.decomposition
    sm = prop.site_model
    @assert length(sequence) == sm.sequence_length "Length of sequence does not match expected SiteModel sequence length (got $(length(sequence.value)), expected $(sm.sequence_length))."

    update_transition_weights!(prop.transition_weights, Δt, sm.μ, dec.λ, dec.V, dec.V⁻¹)

    updated_sequence = Vector{UInt8}(undef, sm.sequence_length)
    for gamma_category in eachindex(sm.μ)
        cumulative_probabilities = cumsum(prop.transition_weights[gamma_category], dims=1)
        for site in sm.variable_sites[gamma_category]
            r = rand(rng)
            old_nucleotide = sequence[site]
            updated_nucleotide = 1
            while r ≥ cumulative_probabilities[updated_nucleotide, old_nucleotide] && updated_nucleotide < N
                updated_nucleotide += 1
            end
            updated_sequence[site] = updated_nucleotide
        end
    end
    return updated_sequence
end