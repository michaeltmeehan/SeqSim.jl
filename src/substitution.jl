
abstract type SubstitutionModel end


@with_kw mutable struct JC <: SubstitutionModel
    frequencies::Vector{Float64}=fill(0.25, 4)
end


@with_kw mutable struct F81 <: SubstitutionModel
    frequencies::Vector{Float64}
end


@with_kw mutable struct HKY <: SubstitutionModel
    tt_ratio::Float64
    frequencies::Vector{Float64}
end


abstract type ClockModel end


mutable struct StrictClock <: ClockModel
    clock_rate::Float64
end


mutable struct SiteModel
    mutation_rate::Float64
    gamma_category_count::Int
    gamma_shape::Float64
    proportion_invariant::Float64
    substitution_model::SubstitutionModel
end


function rate_matrix(π::Vector{Float64}, Q_bare::Matrix{Float64})::Matrix{Float64}
    Q = π .* Q_bare
    for i in 1:4
        Q[i,i] = -sum(Q[:,i])
    end
    return Q
end


function rate_matrix(substitution_model::JC)::Matrix{Float64}
    π = fill(0.25, 4)
    Q_bare = ones(4,4)
    return rate_matrix(π, Q_bare)
end


function rate_matrix(substitution_model::F81)::Matrix{Float64}
    π = substitution_model.frequencies
    Q_bare = ones(4,4)
    return rate_matrix(π, Q_bare)
end


function rate_matrix(substitution_model::HKY)::Matrix{Float64}
    π = substitution_model.frequencies
    κ = substitution_model.tt_ratio
    Q_bare = [0. 1. κ 1.; 1. 0. 1. κ; κ 1. 0. 1.; 1. κ 1. 0.]
    return rate_matrix(π, Q_bare)
end