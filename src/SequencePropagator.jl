function update_transition_weights!(transition_weights::Vector{Matrix{Float64}}, 
                                    Δt::Float64, 
                                    μ::Vector{Float64}, 
                                    λ::SVector{4, Float64}, 
                                    V::SMatrix{4,4,Float64}, 
                                    V⁻¹::SMatrix{4,4,Float64};
                                    validate_weights::Bool=false)::Vector{Matrix{Float64}}
    isfinite(Δt) && Δt >= 0.0 || throw(ArgumentError("Elapsed time Δt must be finite and non-negative."))
    length(transition_weights) == length(μ) || throw(ArgumentError("Number of transition-weight matrices must match number of site-rate categories."))
    for gamma_category in eachindex(transition_weights)
        size(transition_weights[gamma_category]) == (4, 4) || throw(ArgumentError("Transition-weight matrices must be 4x4 for DNA states A, C, G, T."))
        isfinite(μ[gamma_category]) && μ[gamma_category] >= 0.0 || throw(ArgumentError("Site-rate category μ[$gamma_category] must be finite and non-negative."))
        transition_weights[gamma_category] .= V * Diagonal(exp.((μ[gamma_category] * Δt) .* λ)) * V⁻¹
        validate_weights && validate_and_normalize_weights(transition_weights[gamma_category])
    end
    return transition_weights
end


function validate_and_normalize_weights(weights::Vector{Matrix{Float64}})::Vector{Matrix{Float64}}
    for gamma_category in eachindex(weights)
        validate_and_normalize_weights(weights[gamma_category])
    end
    return weights
end

function validate_and_normalize_weights(weights::Matrix{Float64})::Matrix{Float64}
    for col_index in axes(weights, 2)
        col_sum = sum(view(weights, :, col_index))
        if !isapprox(col_sum, 1.0; atol=1e-10)
            @warn "Column $col_index does not sum to 1.0 (sum = $col_sum). Renormalizing."
            weights[:, col_index] ./= col_sum
        end
    end
    return weights
end


struct SequencePropagator{N}
    site_model::SiteModel
    decomposition::RateMatrixDecomposition{N}
    transition_weights::Vector{Matrix{Float64}}
end


function SequencePropagator(site_model::SiteModel)
    Q = rate_matrix(site_model.substitution_model)
    N = size(Q, 1)
    decomposition = decompose(Q)
    transition_weights = [zeros(N, N) for _ in eachindex(site_model.μ)]
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
    length(sequence) == sm.sequence_length || throw(ArgumentError("Length of sequence does not match expected SiteModel sequence length (got $(length(sequence)), expected $(sm.sequence_length))."))
    validate_encoded_sequence(sequence)
    isfinite(Δt) && Δt >= 0.0 || throw(ArgumentError("Elapsed time Δt must be finite and non-negative."))

    update_transition_weights!(prop.transition_weights, Δt, sm.μ, dec.λ, dec.V, dec.V⁻¹)

    updated_sequence = copy(sequence)
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
