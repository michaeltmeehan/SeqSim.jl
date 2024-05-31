# TODO: Consider adding LongSequence option for memory savings (cutoff seems to be around 100 characters)
# seq = dna"ACG"
# To construct sequences see: https://github.com/BioJulia/BioSequences.jl/blob/master/docs/src/construction.md

const nucleotides = [DNA_A, DNA_C, DNA_G, DNA_T]
const nucleotide_idx = Dict(zip(nucleotides, 1:4))


"""
    simulate_sequence(n::Int64, frequencies::Vector{<:Number})::Vector{DNA}

Generate a random sequence of nucleotides for the root of a phylogenetic tree, using the provided frequencies for each nucleotide.
This function leverages `SamplerWeighted` from the BioSequences package, which expects the first three elements of the frequency vector;
the last frequency is automatically calculated to make the total sum to 1.

# Arguments
- `n::Int64`: Length of the sequence to generate.
- `frequencies::Vector{<:Number}`: Frequencies of each nucleotide, corresponding to ['A', 'C', 'G', 'T'].

# Returns
- `Vector{DNA}`: A vector of DNA nucleotides.

# Example
```julia
simulate_sequence(10, [0.1, 0.2, 0.3, 0.4])  # The function automatically calculates the frequency for the fourth nucleotide.
```
"""
function simulate_sequence(n::Int64, frequencies::Vector{<:Number})::Vector{DNA}
    n ≤ 0 && throw(ArgumentError("Sequence length must be a positive integer. Received: $n"))
    length(frequencies) != 4 && throw(ArgumentError("Require four frequencies. Provided $(length(frequencies))"))
    any(frequencies .< 0.) && throw(ArgumentError("Frequencies cannot be negative. Received: $(frequencies)"))
    isapprox(sum(frequencies), 1.0; atol=1e-5) || throw(ArgumentError("Frequencies must sum to 1. Received: $(frequencies)"))
    return rand(SamplerWeighted(nucleotides, frequencies[1:3]), n)  # SamplerWeighted only used first n-1 elements of frequencies (final value is calculated from normalized sum)
end


"""
    simulate_sequences!(tree::RootedTree, seq_length::Int64, site_model::SiteModel)

Simulate DNA sequences for all nodes in a phylogenetic tree based on a specified site model. This function starts
by assigning a sequence to the root of the tree and then propagates these sequences to all descendant nodes,
modifying them according to evolutionary parameters defined by the `SiteModel`.

# Arguments
- `tree::RootedTree`: A phylogenetic tree structure containing nodes linked by branches. Each node represents a taxonomic unit (e.g., species, individual), and branches represent evolutionary relationships and distances.
- `seq_length::Int64`: The length of the DNA sequences to be generated. This length should be a positive integer, indicating the number of nucleotides in each DNA sequence.
- `site_model::SiteModel`: A model encapsulating all evolutionary parameters necessary for sequence simulation. This includes the mutation rate, the distribution of mutation rates across sites (possibly modeled by a gamma distribution), the proportion of invariant sites, and the substitution model used for base changes.

# Details
- The function initializes by checking the validity of the input parameters, including the sequence length and the non-emptiness of the tree.
- It computes the rate matrix from the substitution model provided in the `SiteModel` and decomposes this matrix to use in simulating evolutionary changes along the tree.
- For each node in the tree (starting from the root and excluding it in subsequent iterations), the function simulates the sequence evolution from its parent node using the precomputed evolutionary parameters (mutation rates and transition probabilities).
- Sequences are evolved based on branch lengths, which represent the evolutionary distances between nodes.

# Returns
- The function modifies the `tree` in-place by assigning simulated DNA sequences to each node's `data` dictionary under the key `"sequence"`.

# Examples
```julia
using Phylo  # Assume necessary packages are loaded
# Create a random phylogenetic tree
tree = rand(Nonultrametric(10))  # Create a random tree with 10 nodes
# Define a substitution model
sub_model = HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4])
# Create a site model with no invariant sites and no gamma variation
site_model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0.1, substitution_model=sub_model)
# Simulate sequences
simulate_sequences!(tree, 100, site_model)
```
# Notes

Ensure that the tree has at least one node and that each node, starting from the root, has valid connections to form a well-defined tree structure.
The function throws an ArgumentError if the input conditions are not met, including checks for positive sequence length and non-empty tree structure.

"""
function simulate_sequences!(tree::RootedTree, seq_length::Int64, site_model::SiteModel)
    seq_length ≤ 0 && throw(ArgumentError("Sequence length must be a positive integer. Received: $seq_length"))
    nnodes(tree) == 0 && throw(ArgumentError("The provided tree is empty."))

    @unpack mutation_rate, gamma_category_count, gamma_shape, proportion_invariant, substitution_model = site_model
    π = substitution_model.π
    Q = rate_matrix(substitution_model)
    λ, V, V⁻¹ = decompose(Q)
    P = zeros(4, 4)     # Initialize transition matrix
    μ, weights = assign_rates(seq_length, site_model)

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


@forward Phylogeny.tree simulate_sequences!


"""
    assign_rates(n::Int64, model::SiteModel) -> Tuple

Calculate mutation rates (`μ`) and assign transition weight dictionaries (`weights`) for `n` sites based on the `SiteModel`.
This function adapts to the configuration of the model, handling both invariant and variable rate categories as defined by the gamma distribution.

# Arguments
- `n`: Number of sites.
- `model`: A `SiteModel` encapsulating all necessary evolutionary parameters including rate, shape, category count, and invariance proportion.

# Returns
- `μ`: Either a single value or an array of mutation rates, depending on the site variability.
- `weights`: A dictionary or dictionary of dictionaries containing transition probabilities for each nucleotide.

# Examples
```julia
model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0.1, substitution_model=JC())
rates, weights = assign_rates(100, model)
```
"""
function assign_rates(n::Int64, model::SiteModel)
    @unpack mutation_rate, gamma_category_count, gamma_shape, proportion_invariant, substitution_model = model
    mutation_rate <= 0. && error("Mutation rate must be positive.")
    gamma_shape < 0 && error("Gamma category count must be non-negative.")
    gamma_shape < 0. && error("Gamma shape must be non-negative.")
    (proportion_invariant < 0. || proportion_invariant ≥ 1.) && error("Proportion invariant must be between 0 and 1.")
    if gamma_category_count == 0
         # Handle the case with no gamma categories: assign uniform transition probabilities
        weights = Dict([nucleotides[i] => fill(0.25, 4) for i in 1:4])
        if proportion_invariant == 0.
            μ = mutation_rate
        else
            # Compute mutation rates considering the proportion of invariant sites
            μ = fill(0., n)
            variable_sites = rand(n) .> proportion_invariant
            μ[variable_sites] .= mutation_rate / (1. - proportion_invariant)
        end
    else
        # Compute rates using a gamma distribution for variable sites
        d = Gamma(gamma_shape, mutation_rate / (gamma_shape * (1. - proportion_invariant)))
        q = quantile(d, range(start=1. /(2. * gamma_category_count), step=1. / gamma_category_count, length = gamma_category_count))
        pushfirst!(q, 0.)
        rate_categories = assign_rate_categories(model, n)
        μ = [(rc, q[rc]) for rc in rate_categories]
        weights = Dict([i => Dict([nucleotides[j] => fill(0.25, 4) for j in 1:4]) for i in 2:(1+gamma_category_count)])
    end    
    return μ, weights
end


"""
    propagate_sequence(seq_in::Vector{DNA}, μ::Float64, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}, P::Matrix{Float64}, weights::Dict{DNA, Vector{Float64}})::Vector{DNA} where T <: Number

Simulate the evolution of a DNA sequence over a specified time interval using a constant mutation rate.

# Arguments
- `seq_in`: Input sequence, a vector of `DNA` types representing the genetic sequence.
- `μ`: Constant mutation rate applied to all sites in the sequence.
- `Δt`: Time interval over which the sequence evolves.
- `λ`: Eigenvalues of the rate matrix, used in computing the transition probabilities.
- `V`: Eigenvectors of the rate matrix.
- `V⁻¹`: Inverse of the matrix of eigenvectors.
- `P`: Preallocated matrix for storing computed transition probabilities.
- `weights`: Dictionary mapping each nucleotide to a vector of transition probabilities.

# Returns
- `Vector{DNA}`: The evolved sequence as a vector of `DNA` types.

# Example
```julia
model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0.1, substitution_model=HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4]))
Q = rate_matrix(model)
seq = [DNA_A, DNA_C, DNA_G, DNA_T]
λ, V, V⁻¹ = decompose(Q)
P = zeros(4, 4)
weights = Dict(DNA_A => zeros(4), DNA_C => zeros(4), DNA_G => zeros(4), DNA_T => zeros(4))
μ = 0.1
Δt = 1.0
evolved_seq = propagate_sequence(seq, μ, Δt, λ, V, V⁻¹, P, weights)
```
"""
function propagate_sequence(seq_in::Vector{DNA}, μ::Float64, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}, P::Matrix{Float64}, weights::Dict{DNA, Vector{Float64}})::Vector{DNA} where T <: Number
    seq_out = copy(seq_in)
    compute_transition_weights!(P, weights, μ, Δt, λ, V, V⁻¹)
    for i in eachindex(seq_in)
        seq_out[i] = update_site(seq_in[i], weights)
    end
    return seq_out
end


"""
    propagate_sequence(seq_in::Vector{DNA}, μ::Vector{Float64}, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}, P::Matrix{Float64}, weights::Dict{DNA, Vector{Float64}})::Vector{DNA} where T <: Number

Simulate the evolution of a DNA sequence where a proportion of sites are invariant.

# Arguments
- `seq_in`: Input sequence, a vector of `DNA` types.
- `μ`: Vector of mutation rates for each site in the sequence.
- `Δt`: Time interval over which the sequence evolves.
- `λ`, `V`, `V⁻¹`, `P`, `weights`: Same as in the constant mutation rate method.

# Returns
- `Vector{DNA}`: The evolved sequence.

# Example
```julia
model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, gamma_category_count=4, proportion_invariant=0.1, substitution_model=HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4]))
Q = rate_matrix(model)
seq = [DNA_A, DNA_C, DNA_G, DNA_T]
λ, V, V⁻¹ = decompose(Q)
P = zeros(4, 4)
weights = Dict(DNA_A => zeros(4), DNA_C => zeros(4), DNA_G => zeros(4), DNA_T => zeros(4))
μ = [0.1, 0.2, 0.05, 0.3]  # Different mutation rates for each site
evolved_seq = propagate_sequence(seq, μ, Δt, λ, V, V⁻¹, P, weights)
```
"""
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


"""
    propagate_sequence(seq_in::Vector{DNA}, μ::Vector{Tuple{Int64,Float64}}, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}, P::Matrix{Float64}, weights::Dict{Int64,Dict{DNA, Vector{Float64}}})::Vector{DNA} where T <: Number

Simulate the evolution of a DNA sequence where each site may have a different mutation rate categorized by rate categories. This method allows for complex models where different groups of sites evolve under different evolutionary pressures or mutation rates.

# Arguments
- `seq_in`: Input sequence, a vector of `DNA` types representing the genetic sequence.
- `μ`: Vector of tuples, where each tuple contains a category index and its corresponding mutation rate.
- `Δt`: Time interval over which the sequence evolves.
- `λ`: Eigenvalues of the rate matrix, used in computing the transition probabilities.
- `V`: Eigenvectors of the rate matrix.
- `V⁻¹`: Inverse of the matrix of eigenvectors.
- `P`: Preallocated matrix for storing computed transition probabilities.
- `weights`: Dictionary mapping each category to a dictionary that maps each nucleotide to a vector of transition probabilities.

# Returns
- `Vector{DNA}`: The evolved sequence as a vector of `DNA` types.
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


"""
    compute_transition_weights!(P::Matrix{Float64}, weights::Dict{DNA, Vector{Float64}}, μ::Float64, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}) where T <: Number

Update the transition probability matrix `P` in place based on the provided eigenvalues `λ`, eigenvectors `V`, and their inverse `V⁻¹`. 
The function calculates the matrix exponential of the rate matrix scaled by mutation rate `μ` and time interval `Δt`, and uses this to populate `P`.

This function also populates a dictionary `weights` that maps each nucleotide to its respective transition probabilities. 
    These probabilities are normalized to sum to one for each nucleotide.

# Parameters
- `P`: Preallocated matrix to be filled with the computed transition probabilities.
- `weights`: Dictionary mapping each nucleotide type to a vector of transition probabilities.
- `μ`: Mutation rate per unit time.
- `Δt`: Time interval over which the evolution is considered.
- `λ`: Vector containing the eigenvalues of the rate matrix.
- `V`: Matrix of eigenvectors of the rate matrix.
- `V⁻¹`: Inverse of the matrix of eigenvectors.

# Returns
Nothing. The function modifies `P` and `weights` in place.

# Examples
```julia
P = zeros(4, 4)
weights = Dict(DNA_A => zeros(4), DNA_C => zeros(4), DNA_G => zeros(4), DNA_T => zeros(4))
μ = 0.1
Δt = 1.0
λ = [0.1, -0.1, 0.2, -0.2]
V = rand(4, 4)
V⁻¹ = inv(V)
compute_transition_weights!(P, weights, μ, Δt, λ, V, V⁻¹)
```
"""
@inline function compute_transition_weights!(P::Matrix{Float64}, weights::Dict{DNA, Vector{Float64}}, μ::Float64, Δt::Float64, λ::Vector{T}, V::Matrix{T}, V⁻¹::Matrix{T}) where T <: Number
    # Calculate the matrix exponential of the scaled rate matrix
    P .= V * Diagonal(exp.(μ * λ * Δt)) * V⁻¹

    # Loop over each nucleotide to compute and normalize transition probabilities
    for i in 1:4
        prob = 0.   # Initialize cumulative probability for normalization
        for j in 0:3
            # Calculate transition probability for nucleotide `i` transitioning to `mod_wrap(i+j, 4)`
            prob += P[mod_wrap(i+j, 4), i]
            weights[nucleotides[i]][j+1] = P[mod_wrap(i+j, 4), i]
        end
        # Normalize the probabilities so that they sum to 1
        weights[nucleotides[i]] ./= prob
    end
end


"""
    update_site(nucl::DNA, weights::Dict{DNA, Vector{Float64}}) -> DNA

Update a nucleotide based on given transition weights.

# Arguments
- `nucl`: The current nucleotide.
- `weights`: A dictionary mapping each nucleotide to a vector of transition probabilities.

# Returns
- The updated nucleotide after applying the transition probabilities.

# Example
```julia
nucl = 'A'
weights = Dict('A' => [0.1, 0.2, 0.3, 0.4])
updated_nucl = update_site(nucl, weights)
```
"""
@inline function update_site(nucl::DNA, weights::Dict{DNA, Vector{Float64}})
    !haskey(weights, nucl) && throw(ArgumentError("No transition weights available for nucleotide '$nucl'"))

    # Generate a random number for selecting the transition
    r = rand()
    prob = weights[nucl][1]

    # Check if the nucleotide remains unchanged
    if r ≤ prob
        return nucl
    else
        # Find the new nucleotide state based on cumulative weights
        nucl_idx = nucleotide_idx[nucl]
        idx = 1
        while r ≥ prob
            idx += 1
            prob += weights[nucl][idx]
        end
        idx ≤ 1 && throw(RuntimeError("Failed to find a valid transition for nucleotide '$nucl'. This might indicate an issue with the weights."))
    end
    # Return the updated nucleotide state, adjusting for index wrapping
    return nucleotides[mod_wrap(idx + nucl_idx - 1, 4)]
end