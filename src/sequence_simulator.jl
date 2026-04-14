
function rand_seq_int(rng::AbstractRNG, sequence_length::Int, frequencies::Vector{Float64})
    sequence_length > 0 || throw(ArgumentError("Sequence length must be positive."))
    validate_frequencies(frequencies)
    d = Categorical(frequencies)
    seq = rand(rng, d, sequence_length)
    return UInt8.(seq)  # Convert to UInt8
end


function rand_seq_int(rng::AbstractRNG, model::SiteModel)
    π = get_frequencies(model.substitution_model)
    return rand_seq_int(rng, model.sequence_length, π)
end


function rand_seq_int(model::SiteModel)
    return rand_seq_int(Random.GLOBAL_RNG, model)
end


function rand_seq(rng::AbstractRNG, sequence_length::Int; frequencies::Vector{Float64} = fill(0.25, 4), taxon=nothing, time=nothing)::Sequence
    return Sequence(rand_seq_int(rng, sequence_length, frequencies), taxon=taxon, time=time)
end


function rand_seq(rng::AbstractRNG, model::SiteModel; taxon=nothing, time=nothing)::Sequence
    return Sequence(rand_seq_int(rng, model), taxon=taxon, time=time)
end


function rand_seq(model::SiteModel; taxon=nothing, time=nothing)::Sequence
    return Sequence(rand_seq_int(model), taxon=taxon, time=time)
end
