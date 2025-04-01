
"""
    update_transition_weights!(transition_weights::Vector{Matrix{Float64}}, 
                               Δt::Float64, 
                               μ::Vector{Float64}, 
                               λ::SVector{4, Float64}, 
                               V::SMatrix{4,4,Float64}, 
                               V⁻¹::SMatrix{4,4,Float64})::Vector{Matrix{Float64}}

Updates the transition weight matrices in-place for a set of categories based on the given parameters.

# Arguments
- `transition_weights::Vector{Matrix{Float64}}`: A vector of matrices representing the transition weights for each category. These matrices are updated in-place.
- `Δt::Float64`: The time interval over which the transition weights are updated.
- `μ::Vector{Float64}`: A vector of mutation rates, one for each category.
- `λ::SVector{4, Float64}`: A static vector representing the eigenvalues of the rate matrix.
- `V::SMatrix{4,4,Float64}`: A static matrix representing the eigenvectors of the rate matrix.
- `V⁻¹::SMatrix{4,4,Float64}`: A static matrix representing the inverse of the eigenvector matrix `V`.

# Returns
- `Vector{Matrix{Float64}}`: The updated `transition_weights` vector (same object as the input, modified in-place).

# Details
For each category, the function computes the updated transition weight matrix using the formula:

    `transition_weights[cat] .= V * Diagonal(exp.(μ[cat] * λ * Δt)) * V⁻¹`

This formula applies the matrix exponential of the scaled eigenvalues to compute the transition probabilities over the given time interval `Δt`.

# Notes
- The function modifies the `transition_weights` vector in-place, so the input object is directly updated.
- The use of static arrays (`SVector` and `SMatrix`) for `λ`, `V`, and `V⁻¹` ensures efficient computation for small fixed-size matrices.
"""
function update_transition_weights!(transition_weights::Vector{Matrix{Float64}}, 
                                    Δt::Float64, 
                                    μ::Vector{Float64}, 
                                    λ::SVector{4, Float64}, 
                                    V::SMatrix{4,4,Float64}, 
                                    V⁻¹::SMatrix{4,4,Float64})::Vector{Matrix{Float64}}
    for cat in eachindex(transition_weights)
        transition_weights[cat] = V * Diagonal(exp.((μ[cat] * Δt) .* λ)) * V⁻¹    # TODO: Replace with V * (exp.(μ[cat] * λ * Δt) .* V⁻¹)
        # if !all(sum(transition_weights[cat], dims=1) .≈ 1.0)
        #     @warn "Transition weights for category $cat do not sum to 1.0."
        #     for col in eachcol(transition_weights[cat])
        #         col_sum = sum(col)
        #         if !isapprox(col_sum, 1.0; atol=1e-10)
        #             @warn "Column $col in transition_weights[$cat] does not sum to 1.0 (sum = $col_sum). Renormalizing."
        #             transition_weights[cat][:, col] ./= col_sum
        #         end
        #     end
        # end
    end
    return transition_weights
end


"""
    update_site(nucl_in::UInt8, weights::Matrix{Float64})::UInt8

Updates a nucleotide site based on a given transition probability matrix.

# Arguments
- `nucl_in::UInt8`: The current nucleotide represented as an unsigned 8-bit integer.
- `weights::Matrix{Float64}`: A matrix of transition probabilities where each entry `weights[i, j]` 
  represents the probability of transitioning from nucleotide `j` to nucleotide `i`.

# Returns
- `UInt8`: The updated nucleotide after applying the transition probabilities.

# Details
The function generates a random number `r` and iteratively accumulates probabilities from the 
`weights` matrix until `r` is less than the accumulated probability. The index of the accumulated 
probability determines the new nucleotide.

# Notes
- The input `weights` matrix is assumed to be properly normalized such that the sum of probabilities 
  for each column equals 1.
- The function uses 1-based indexing, consistent with Julia's array indexing.
"""
@inline function update_site(nucl_in::UInt8, weights::Matrix{Float64})::UInt8
    r = rand()
    prob = weights[1, nucl_in]
    nucl_out = UInt8(1)
    while r ≥ prob
        nucl_out += 1
        prob += weights[nucl_out, nucl_in]
    end
    return nucl_out
end


"""
    update_sequence!(transition_weights::Vector{Matrix{Float64}}, 
                     seq_in::Vector{UInt8}, 
                     Δt::Float64, 
                     μ::Vector{Float64},
                     variable_sites::Vector{Vector{Int64}},
                     λ::SVector{4, Float64},
                     V::SMatrix{4,4,Float64},
                     V⁻¹::SMatrix{4,4,Float64})::Vector{UInt8}

Updates a nucleotide sequence based on a given set of transition weights and evolutionary parameters.

# Arguments
- `transition_weights::Vector{Matrix{Float64}}`: A vector of transition weight matrices, one for each category.
- `seq_in::Vector{UInt8}`: The input nucleotide sequence represented as a vector of `UInt8` values.
- `Δt::Float64`: The time interval over which the sequence evolves.
- `μ::Vector{Float64}`: A vector of mutation rates for each category.
- `variable_sites::Vector{Vector{Int64}}`: A vector where each element is a list of indices representing variable sites for a specific category.
- `λ::SVector{4, Float64}`: A static vector of eigenvalues for the substitution model.
- `V::SMatrix{4,4,Float64}`: A static matrix of eigenvectors for the substitution model.
- `V⁻¹::SMatrix{4,4,Float64}`: The inverse of the eigenvector matrix `V`.

# Returns
- `Vector{UInt8}`: The updated nucleotide sequence after applying the transition weights and evolutionary parameters.

# Details
This function simulates the evolution of a nucleotide sequence over a given time interval (`Δt`). It updates the sequence by iterating over each category of sites and applying the corresponding transition weights to the variable sites. The transition weights are updated using the `update_transition_weights!` function, which incorporates the evolutionary parameters (`λ`, `V`, `V⁻¹`, and `μ`).

Each site in the sequence is updated using the `update_site` function, which determines the new nucleotide state based on the current state and the transition probabilities.

# Note
The function assumes that the input sequence (`seq_in`) and the transition weights are properly initialized and consistent with the provided evolutionary parameters.
"""
function update_sequence!(transition_weights::Vector{Matrix{Float64}}, 
                          seq_in::Vector{UInt8}, 
                          Δt::Float64, 
                          μ::Vector{Float64},
                          variable_sites::Vector{Vector{Int64}},
                          λ::SVector{4, Float64},
                          V::SMatrix{4,4,Float64},
                          V⁻¹::SMatrix{4,4,Float64})::Vector{UInt8}
    seq_out = copy(seq_in)
    transition_weights = update_transition_weights!(transition_weights, Δt, μ, λ, V, V⁻¹)
    for cat in eachindex(μ)
        weights = transition_weights[cat]
        for idx in variable_sites[cat]
            seq_out[idx] = update_site(seq_in[idx], weights)
        end
    end
    return seq_out
end