
const nucleotides = ['a', 'c', 'g', 't']
const nucleotide_idx = Dict(zip(nucleotides, 1:4))


"""
    simulate_sequence(n::Int64, frequencies::Vector{Float64})::Vector{Char}

Generate a random sequence of nucleotides for the root of a phylogenetic tree.

# Arguments
- `n::Int64`: Length of the sequence to generate.
- `frequencies::Vector{Float64}`: Frequencies of each nucleotide, corresponding to ['a', 'c', 'g', 't'].

# Returns
- `Vector{Char}`: A vector of characters representing the sequence of nucleotides.

# Example
```julia
simulate_sequence(10, [0.1, 0.2, 0.3, 0.4])
"""
function simulate_sequence(n::Int64, frequencies::Vector{Float64})::Vector{Char}
    if sum(frequencies) ≈ 1.0
        return sample(nucleotides, Weights(frequencies), n, replace=true)
    else
        throw(ArgumentError("Frequencies must sum to 1"))
    end
end


function simulate_sequences!(tree::RootedTree, seq_length::Int64, site_model::SiteModel)
    @unpack mutation_rate, gamma_category_count, gamma_shape, substitution_model = site_model
    μ = mutation_rate
    π = substitution_model.frequencies
    Q = rate_matrix(substitution_model)
    D, V, V⁻¹ = decompose(Q)
    getroot(tree).data["sequence"] = simulate_sequence(seq_length, π)

    for node in Iterators.drop(traversal(tree, preorder), 1)
        parent = getparent(tree, node)
        Δt = getbranch(tree, parent, node).length
        node.data["sequence"] = propagate_sequence(parent.data["sequence"], μ, Δt, D, V, V⁻¹)
    end
end


function compute_transition_weights(μ::Float64, Δt::Float64, D::Matrix{T}, V::Matrix{T}, V⁻¹::Matrix{T}) where T <: Number
    P = V * exp(μ * D * Δt) * V⁻¹
    maximum(imag.(P)) > 1e-10 && @warn "Transition matrix has imaginary parts > 1e-10"
    P = real(P)
    weights = Dict{Char, Vector}()
    for i in 1:4
        weights[nucleotides[i]] = cumsum([P[mod2(i+j, 4), i] for j in 0:3]) # Returned as cyclic permutations
        weights[nucleotides[i]] ./= weights[nucleotides[i]][end]
    end
    return weights
end


function propagate_sequence(seq_in::Vector{Char}, μ::Float64, Δt::Float64, D::Matrix{T}, V::Matrix{T}, V⁻¹::Matrix{T})::Vector{Char} where T <: Number
    seq_out = copy(seq_in)
    weights = compute_transition_weights(μ, Δt, D, V, V⁻¹)
    for i in eachindex(seq_in)
        nucl = seq_in[i]
        r = rand()
        if r > weights[nucl][1]
            nucl_idx = nucleotide_idx[nucl]
            idx = findfirst(x -> x >= r, weights[nucl])
            seq_out[i] = nucleotides[mod2(idx + nucl_idx - 1, 4)]
        end
    end
    return seq_out
end