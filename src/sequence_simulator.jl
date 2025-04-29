
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