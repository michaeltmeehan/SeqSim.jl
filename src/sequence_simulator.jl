
function rand_seq_int(rng::AbstractRNG, sequence_length::Int; frequencies::Vector{Float64} = fill(0.25, 4))
    @assert length(frequencies) == 4 "Frequencies must be a vector of length 4."
    @assert sum(frequencies) ≈ 1.0 "Frequencies must sum to 1. Received sum = $(sum(frequencies))"
    @assert all(frequencies .>= 0) "Frequencies must be non-negative. Received frequencies = $frequencies"
    sequence_length > 0 || throw(ArgumentError("Sequence length must be positive."))
    d = Categorical(frequencies)
    seq = rand(rng, d, sequence_length)
    return UInt8.(seq)  # Convert to UInt8
end


function rand_seq(rng::AbstractRNG, sequence_length::Int; frequencies::Vector{Float64} = fill(0.25, 4), taxon=nothing, time=nothing)::Sequence
    return Sequence(rand_seq_int(rng, sequence_length; frequencies=frequencies), taxon=taxon, time=time)
end


function rand_seq(sequence_length::Int; frequencies::Vector{Float64}=fill(0.25, 4), taxon=nothing, time=nothing)::Sequence
    return Sequence(rand_seq_int(Random.GLOBAL_RNG, sequence_length; frequencies=frequencies), taxon=taxon, time=time)
end


function update_transition_weights!(transition_weights::Vector{Matrix{Float64}}, 
                                    Δt::Float64, 
                                    μ::Vector{Float64}, 
                                    λ::SVector{4, Float64}, 
                                    V::SMatrix{4,4,Float64}, 
                                    V⁻¹::SMatrix{4,4,Float64};
                                    validate_weights::Bool=false)::Vector{Matrix{Float64}}
    for gamma_category in eachindex(transition_weights)
        transition_weights[gamma_category] .= V * Diagonal(exp.((μ[gamma_category] * Δt) .* λ)) * V⁻¹
        validate_weights && validate_and_normalize_weights(transition_weights[gamma_category])
    end
    return transition_weights
end


function validate_and_normalize_weights(weights::Vector{Matrix{Float64}})::Vector{Matrix{Float64}}
    for gamma_category in eachindex(weights)
        for col in eachcol(weights[gamma_category])
            col_sum = sum(col)
            if !isapprox(col_sum, 1.0; atol=1e-10)
                @warn "Column $col in transition_weights[$gamma_category] does not sum to 1.0 (sum = $col_sum). Renormalizing."
                weights[gamma_category][:, col] ./= col_sum
            end
        end
    end
    return weights
end


function update_sequence!(rng::AbstractRNG,
                          sequence::Vector{UInt8}, 
                          transition_weights::Vector{Matrix{Float64}}, 
                          Δt::Float64, 
                          μ::Vector{Float64},
                          variable_sites::Vector{Vector{Int64}},
                          λ::SVector{4, Float64},
                          V::SMatrix{4,4,Float64},
                          V⁻¹::SMatrix{4,4,Float64})::Vector{UInt8}
    update_transition_weights!(transition_weights, Δt, μ, λ, V, V⁻¹)
    for gamma_category in eachindex(μ)
        cumulative_probabilities = cumsum(transition_weights[gamma_category], dims=1)
        @inbounds for site in variable_sites[gamma_category]
            r = rand(rng)
            old_nucleotide = sequence[site]
            updated_nucleotide = 1
            while r ≥ cumulative_probabilities[updated_nucleotide, old_nucleotide]
                updated_nucleotide += 1
            end
            sequence[site] = updated_nucleotide
        end
    end
    return sequence
end

struct SequenceSimulator{N}
    site_model::SiteModel
    decomposition::RateMatrixDecomposition{N}
    transition_weights::Vector{Matrix{Float64}}
end


function SequenceSimulator(site_model::SiteModel)
    Q = rate_matrix(site_model.substitution_model)
    N = size(Q, 1)
    decomposition = decompose(Q)
    transition_weights = [zeros(N, N) for _ in 1:site_model.gamma_category_count]
    return SequenceSimulator{N}(site_model, decomposition, transition_weights)
end


# function (sim::SequenceSimulator{N})(rng::AbstractRNG, sequence::Sequence, Δt::Float64)
#     sm = sim.site_model
#     dec = sim.decomposition

#     update_transition_weights!(sim.transition_weights, Δt, sm.μ, dec.λ, dec.V, dec.V⁻¹)

#     for gamma_category in eachindex(sm.μ)
#         cumulative_probabilities = cumsum(sim.transition_weights[gamma_category], dims=1)
#         @inbounds for site in sm.variable_sites[gamma_category]
#             r = rand(rng)
#             old_nucleotide = sequence.value[site]
#             updated_nucleotide = 1
#             while r ≥ cumulative_probabilities[updated_nucleotide, old_nucleotide]
#                 updated_nucleotide += 1
#             end
#             sequence.value[site] = updated_nucleotide
#         end
#     end

#     if !isnothing(sequence.time)
#         sequence.time += Δt
#     end

#     return sequence
# end



function update_sequence!(rng::AbstractRNG, 
                          sequence::Sequence, 
                          site_model::SiteModel, 
                          Δt::Float64)

    # Extract parameters from the site model
    variable_sites = site_model.variable_sites
    μ = site_model.μ
    λ = site_model.substitution_model.λ
    V = site_model.substitution_model.V
    V⁻¹ = site_model.substitution_model.V⁻¹

    # Initialize transition weights for each category
    transition_weights = [zeros(4, 4) for _ in 1:site_model.gamma_category_count]

    # Update the sequence using the transition weights and parameters
    updated_sequence = update_sequence!(rng, sequence.value, transition_weights, Δt, μ, variable_sites, λ, V, V⁻¹)

    # Update sequence time if available
    if !isnothing(sequence.time)
        sequence.time += Δt
    end

    return updated_sequence
end