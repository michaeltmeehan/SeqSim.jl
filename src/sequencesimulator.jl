
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
    n ≤ 0 && throw(ArgumentError("Sequence length must be a positive integer. Received: $n"))
    length(frequencies) != 4 && throw(ArgumentError("Require four frequencies. Provided $(length(frequencies))"))
    any(frequencies .< 0.) && throw(ArgumentError("Frequencies cannot be negative. Received: $(frequencies)"))
    isapprox(sum(frequencies), 1.0; atol=1e-5) || throw(ArgumentError("Frequencies must sum to 1. Received: $(frequencies)"))
    return sample(nucleotides, Weights(frequencies), n, replace=true)
end


"""
    simulate_sequences!(tree::RootedTree, seq_length::Int64, site_model::SiteModel)

Simulate DNA sequences for each node in a phylogenetic tree based on a given site model. 
This function first simulates a sequence for the root based on the equilibrium frequencies and then
propagates this sequence to all descendants, accounting for evolutionary changes dictated by the site model.

# Arguments
- `tree::RootedTree`: The phylogenetic tree structure.
- `seq_length::Int64`: Length of the sequence to be generated for each node.
- `site_model::SiteModel`: A model encapsulating the mutation rate, substitution model, and other evolutionary parameters.

# Example
```julia
using Phylo
tree = rand(Nonultrametric(50)) # Create a random tree using the Phylo package
site_model = SiteModel(
    mutation_rate = 1.0,
    gamma_category_count = 4,
    gamma_shape = 0.5,
    proportion_invariant = 0.1,
    substitution_model = JC()
)
simulate_sequences!(tree, 100, site_model)
"""
function simulate_sequences!(tree::RootedTree, seq_length::Int64, site_model::SiteModel)
    seq_length ≤ 0 && throw(ArgumentError("Sequence length must be a positive integer. Received: $seq_length"))
    nnodes(tree) == 0 && throw(ArgumentError("The provided tree is empty."))

    @unpack mutation_rate, gamma_category_count, gamma_shape, substitution_model = site_model
    μ = mutation_rate
    π = substitution_model.π
    Q = rate_matrix(substitution_model)
    λ, V, V⁻¹ = decompose(Q)
    P = zeros(4, 4)     # Initialize transition matrix
    root = getroot(tree)
    isnothing(root) && throw(ArgumentError("The tree has no root."))
    root.data["sequence"] = simulate_sequence(seq_length, π)

    for node in Iterators.drop(traversal(tree, preorder), 1)    # Skipping the root node
        parent = getparent(tree, node)
        isnothing(parent) && throw(ArgumentError("A node has no parent, which is unexpected in a rooted tree."))
        Δt = getbranch(tree, parent, node).length
        node.data["sequence"] = propagate_sequence(parent.data["sequence"], μ, Δt, λ, V, V⁻¹, P)
    end
end


"""
    compute_transition_weights!(P::Matrix{Float64}, μ::Float64, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}) where T <: Number

Mutates the transition probability matrix `P` in place based on the rate matrix's eigen-decomposition and evolutionary parameters.

# Arguments
- `P`: The transition probability matrix to be computed in place.
- `μ`: Mutation rate per unit time.
- `Δt`: Time duration over which the evolution occurs (branch length in phylogenetic terms).
- `λ`: Eigenvalues of the rate matrix.
- `V`: Matrix of eigenvectors of the rate matrix.
- `V⁻¹`: Inverse of the matrix of eigenvectors of the rate matrix.

# Returns
- `weights`: A dictionary mapping each nucleotide to its cumulative transition probabilities.

# Notes
- This function assumes that the eigen-decomposition (`D`, `V`, `V⁻¹`) is precomputed and passed as parameters.
- Warnings are issued if the computed transition matrix contains significant imaginary parts, suggesting numerical issues.

# Example
```julia
# Precomputed matrices for a given rate matrix Q
λ, V, V⁻¹ = decompose(Q)
P = zeros(Float64, 4, 4)
weights = compute_transition_weights!(P, 1.0, 0.05, λ, V, V⁻¹)
"""
function compute_transition_weights!(P::Matrix{Float64}, μ::Float64, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}) where T <: Number
    P .= V * diagm(exp.(μ * λ * Δt)) * V⁻¹
    weights = Dict{Char, Vector}()
    for i in 1:4
        weights[nucleotides[i]] = cumsum([P[mod_wrap(i+j, 4), i] for j in 0:3]) # Returned as cyclic permutations
        weights[nucleotides[i]] ./= weights[nucleotides[i]][end]
    end
    return weights
end


"""
    propagate_sequence(seq_in::Vector{Char}, μ::Float64, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}) where T <: Number -> Vector{Char}

Propagate a sequence through evolutionary time, applying nucleotide transitions based on precomputed transition probabilities.

# Arguments
- `seq_in`: The input sequence of nucleotides.
- `μ`: Mutation rate per unit time.
- `Δt`: Evolutionary time interval (e.g., branch length in a phylogenetic tree).
- `λ`, `V`, `V⁻¹`: Precomputed matrices from the eigen-decomposition of the rate matrix.

# Returns
- A new sequence with mutations applied based on the computed transition probabilities.

# Example
```julia
seq = ['a', 'c', 'g', 't']
μ = 0.1
Δt = 0.05
λ, V, V⁻¹ = decompose(Q)
new_seq = propagate_sequence(seq, μ, Δt, λ, V, V⁻¹)
"""
function propagate_sequence(seq_in::Vector{Char}, μ::Float64, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}, P::Matrix{Float64})::Vector{Char} where T <: Number
    seq_out = copy(seq_in)
    weights = compute_transition_weights!(P, μ, Δt, λ, V, V⁻¹)
    for i in eachindex(seq_in)
        nucl = seq_in[i]
        r = rand()
        if r > weights[nucl][1]
            nucl_idx = nucleotide_idx[nucl]
            idx = findfirst(x -> x >= r, weights[nucl])
            seq_out[i] = nucleotides[mod_wrap(idx + nucl_idx - 1, 4)]
        end
    end
    return seq_out
end