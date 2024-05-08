
const nucleotides = [DNA_A, DNA_C, DNA_G, DNA_T]
const nucleotide_idx = Dict(zip(nucleotides, 1:4))

const P = [1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 1. 0.; 0. 0. 0. 1.]
# const weights = Dict([nucleotides[i] => P[:, i] for i in 1:4])


"""
    simulate_sequence(n::Int64, frequencies::Vector{Float64})::Vector{Char}

Generate a random sequence of nucleotides for the root of a phylogenetic tree.

# Arguments
- `n::Int64`: Length of the sequence to generate.
- `frequencies::Vector{Float64}`: Frequencies of each nucleotide, corresponding to ['a', 'c', 'g', 't'].

# Returns
- `Vector{DNA}`: A vector of characters representing the sequence of nucleotides.

# Example
```julia
simulate_sequence(10, [0.1, 0.2, 0.3, 0.4])
"""
function simulate_sequence(n::Int64, frequencies::Vector{Float64})::Vector{DNA}
    n ≤ 0 && throw(ArgumentError("Sequence length must be a positive integer. Received: $n"))
    length(frequencies) != 4 && throw(ArgumentError("Require four frequencies. Provided $(length(frequencies))"))
    any(frequencies .< 0.) && throw(ArgumentError("Frequencies cannot be negative. Received: $(frequencies)"))
    isapprox(sum(frequencies), 1.0; atol=1e-5) || throw(ArgumentError("Frequencies must sum to 1. Received: $(frequencies)"))
    return rand(SamplerWeighted(nucleotides, frequencies[1:3]), n)
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

    @unpack mutation_rate, gamma_category_count, gamma_shape, proportion_invariant, substitution_model = site_model
    π = substitution_model.π
    Q = rate_matrix(substitution_model)
    λ, V, V⁻¹ = decompose(Q)
    P = zeros(4, 4)     # Initialize transition matrix

    if gamma_category_count == 0
        weights = Dict([nucleotides[i] => P[:, i] for i in 1:4])
        if proportion_invariant == 0.
            μ = mutation_rate
        else
            μ = fill(0., seq_length)
            variable_sites = rand(seq_length) .> proportion_invariant
            μ[variable_sites] .= mutation_rate / (1. - proportion_invariant)
            println("μ: ", μ)
        end
    else
        d = Gamma(gamma_shape, mutation_rate / (gamma_shape * (1. - proportion_invariant)))
        q = quantile(d, range(start=1. /(2. * gamma_category_count), step=1. / gamma_category_count, length = gamma_category_count))
        pushfirst!(q, 0.)
        rate_categories = rate_cat(site_model, seq_length)
        μ = [(rc, q[rc]) for rc in rate_categories]
        weights = Dict([i => Dict([nucleotides[j] => P[:, j] for j in 1:4]) for i in 2:(1+gamma_category_count)])
    end

    root = getroot(tree)
    isnothing(root) && throw(ArgumentError("The tree has no root."))
    root.data["sequence"] = simulate_sequence(seq_length, π)

    for node in Iterators.drop(traversal(tree, preorder), 1)    # Skipping the root node
        parent = getparent(tree, node)
        isnothing(parent) && throw(ArgumentError("A node has no parent, which is unexpected in a rooted tree."))
        Δt = getbranch(tree, parent, node).length
        node.data["sequence"] = propagate_sequence(parent.data["sequence"], μ, Δt, λ, V, V⁻¹, P, weights)
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

# Example
```julia
# Precomputed matrices for a given rate matrix Q
λ, V, V⁻¹ = decompose(Q)
P = zeros(Float64, 4, 4)
weights = compute_transition_weights!(P, 1.0, 0.05, λ, V, V⁻¹)
"""
function compute_transition_weights!(P::Matrix{Float64}, weights::Dict{DNA, Vector{Float64}}, μ::Float64, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}) where T <: Number
    P .= V * Diagonal(exp.(μ * λ * Δt)) * V⁻¹
    for i in 1:4
        prob = 0.
        for j in 0:3
            prob += P[mod_wrap(i+j, 4), i]
            weights[nucleotides[i]][j+1] = P[mod_wrap(i+j, 4), i]
        end
        weights[nucleotides[i]] ./= prob
    end
end


function compute_transition_weights!(P::Matrix{Float64}, weights::Dict{Int64, Dict{DNA, Vector{Float64}}}, μ::Float64, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}) where T <: Number
    for (_, val) in weights
        compute_transition_weights!(P, val, μ, Δt, λ, V, V⁻¹)
    end 
end


"""
    propagate_sequence(seq_in::Vector{DNA}, μ::Float64, Δt::Float64, proportion_invariant::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}, weights::Dict{DNA, Vector{Float64}}) where T <: Number -> Vector{Char}

Propagate a sequence through evolutionary time, applying nucleotide transitions based on precomputed transition probabilities.

# Arguments
- `seq_in`: The input sequence of nucleotides.
- `μ`: Mutation rate per unit time.
- `Δt`: Evolutionary time interval (e.g., branch length in a phylogenetic tree).
- `proportion_invariant`: Proportion of sites that are invariant
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
function propagate_sequence(seq_in::Vector{DNA}, μ::Float64, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}, P::Matrix{Float64}, weights::Dict{DNA, Vector{Float64}})::Vector{DNA} where T <: Number
    seq_out = copy(seq_in)
    compute_transition_weights!(P, weights, μ, Δt, λ, V, V⁻¹)
    for i in eachindex(seq_in)
        seq_out[i] = update_site(seq_in[i], weights)
    end
    return seq_out
end


function propagate_sequence(seq_in::Vector{DNA}, μ::Vector{Float64}, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}, P::Matrix{Float64}, weights::Dict{DNA, Vector{Float64}})::Vector{DNA} where T <: Number
    seq_out = copy(seq_in)
    compute_transition_weights!(P, weights, maximum(μ), Δt, λ, V, V⁻¹)
    for i in eachindex(seq_in)
        if μ[i] > 0.
            seq_out[i] = update_site(seq_in[i], weights)
        end
    end
    return seq_out
end


function update_site(nucl::DNA, weights::Dict{DNA, Vector{Float64}})
    r = rand()
    prob = weights[nucl][1]
    if r ≤ prob
        return nucl
    else
        nucl_idx = nucleotide_idx[nucl]
        idx = 1
        while r ≥ prob
            idx += 1
            prob += weights[nucl][idx]
        end
    end
    return nucleotides[mod_wrap(idx + nucl_idx - 1, 4)]
end


"""
    propagate_sequence(seq_in::Vector{DNA}, μ::Vector{Float64}, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}, weights::Dict{DNA, Vector{Float64}}) where T <: Number -> Vector{Char}

Propagate a sequence through evolutionary time, applying nucleotide transitions based on precomputed transition probabilities.

# Arguments
- `seq_in`: The input sequence of nucleotides.
- `μ`: Site-specific mutation rate per unit time.
- `Δt`: Evolutionary time interval (e.g., branch length in a phylogenetic tree).
- `λ`, `V`, `V⁻¹`: Precomputed matrices from the eigen-decomposition of the rate matrix.

# Returns
- A new sequence with mutations applied based on the computed transition probabilities.

# Example
```julia
seq = ['a', 'c', 'g', 't']
μ = [0., 0.1, 0., 0.2]
Δt = 0.05
λ, V, V⁻¹ = decompose(Q)
new_seq = propagate_sequence(seq, μ, Δt, λ, V, V⁻¹)
"""
function propagate_sequence(seq_in::Vector{DNA}, μ::Vector{Tuple{Int64,Float64}}, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}, P::Matrix{Float64}, weights::Dict{Int64,Dict{DNA, Vector{Float64}}})::Vector{DNA} where T <: Number
    seq_out = copy(seq_in)
    for (cat, rate) in unique(μ)
        rate > 0. && compute_transition_weights!(P, weights[cat], rate, Δt, λ, V, V⁻¹)
    end
    for i in eachindex(seq_in)
        if μ[i][2] > 0.
            seq_out[i] = update_site(seq_in[i], weights[μ[i][1]])
        end
    end
    return seq_out
end