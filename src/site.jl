"""
    SiteModel(sequence_length::Int64, mutation_rate::Float64, gamma_category_count::Int64, gamma_shape::Float64, proportion_invariant::Float64, substitution_model::SubstitutionModel)

Represents a model for site-specific mutation rates in a sequence.

# Arguments
- `sequence_length::Int64`: The total number of sites in the sequence.
- `mutation_rate::Float64`: The overall mutation rate applied to each site.
- `gamma_category_count::Int64`: The number of discrete categories for the gamma distribution of mutation rates.
- `gamma_shape::Float64`: The shape parameter for the gamma distribution.
- `proportion_invariant::Float64`: The proportion of sites that are invariant (i.e., have a mutation rate of zero).
- `substitution_model::SubstitutionModel`: The substitution model used for mutations.

# Fields
- `variable_sites::Vector{Vector{Int64}}`: A vector where each element is a vector of site indices corresponding to a specific mutation rate category.
- `μ::Vector{Float64}`: A vector of mutation rates corresponding to each category of sites.

# Example
```julia
sub_model = SubstitutionModel(...)  # Define your substitution model
site_model = SiteModel(100, 0.01, 4, 0.5, 0.1, sub_model)
"""
struct SiteModel
    sequence_length::Int64
    mutation_rate::Float64
    gamma_category_count::Int64
    gamma_shape::Float64
    proportion_invariant::Float64
    substitution_model::SubstitutionModel
    variable_sites::Vector{Vector{Int64}}
    μ::Vector{Float64}
end

function SiteModel(
    sequence_length::Int64,
    mutation_rate::Float64,
    gamma_category_count::Int64,
    gamma_shape::Float64,
    proportion_invariant::Float64,
    substitution_model::SubstitutionModel
)
    # Input validation
    sequence_length <= 0 && throw(ArgumentError("Sequence length must be positive."))
    mutation_rate <= 0.0 && throw(ArgumentError("Mutation rate must be positive."))
    gamma_category_count < 0 && throw(ArgumentError("Gamma category count must be non-negative."))
    gamma_shape < 0.0 && throw(ArgumentError("Gamma shape must be non-negative."))
    (proportion_invariant < 0.0 || proportion_invariant >= 1.0) && throw(ArgumentError("Proportion invariant must be between 0 and 1."))

    # Assign rates and variable sites
    variable_sites, μ = assign_rates(sequence_length, proportion_invariant, mutation_rate, gamma_shape, gamma_category_count)

    return SiteModel(
        sequence_length,
        mutation_rate,
        gamma_category_count,
        gamma_shape,
        proportion_invariant,
        substitution_model,
        variable_sites,
        μ
    )
end


"""
    assign_rates(sequence_length::Int, proportion_invariant::Float64, mutation_rate::Float64, gamma_shape::Float64, gamma_category_count::Int64) -> (Vector{Vector{Int}}, Vector{Float64})

Assigns mutation rates to sites in a sequence, accounting for invariant sites and discretized gamma-distributed rate variation.

# Arguments
- `sequence_length::Int`: Total number of sites in the sequence.
- `proportion_invariant::Float64`: Proportion of sites that are invariant (i.e., have zero mutation rate), specified as a fraction between 0 and 1.
- `mutation_rate::Float64`: Baseline mutation rate applied to variable sites.
- `gamma_shape::Float64`: Shape parameter of the gamma distribution used to model rate heterogeneity among sites.
- `gamma_category_count::Int64`: Number of discrete categories used to approximate the gamma distribution.

# Returns
- `variable_sites::Vector{Vector{Int}}`: A vector where each element is a list of site indices corresponding to a specific mutation rate category.
- `μ::Vector{Float64}`: A vector of mutation rates corresponding to each rate category.

# Description
This function partitions the sites of a sequence into invariant and variable categories based on the specified `proportion_invariant`. Variable sites are further divided into discrete categories to approximate a gamma distribution of mutation rates, allowing for modeling of rate heterogeneity among sites.

# Example
```julia
sequence_length = 1000
proportion_invariant = 0.2
mutation_rate = 1e-3
gamma_shape = 0.5
gamma_category_count = 4

variable_sites, μ = assign_rates(sequence_length, proportion_invariant, mutation_rate, gamma_shape, gamma_category_count)

In this example, a sequence of 1000 sites is partitioned such that 20% of the sites are invariant, and the remaining sites are divided into 4 categories with mutation rates drawn from a gamma distribution with a shape parameter of 0.5. 
"""
function assign_rates(
    sequence_length::Int,
    proportion_invariant::Float64,
    mutation_rate::Float64,
    gamma_shape::Float64,
    gamma_category_count::Int64
)
    # Determine the number of invariant and variable sites
    num_invariant_sites = Int(floor(proportion_invariant * sequence_length))
    num_variable_sites = sequence_length - num_invariant_sites

    # Initialize variable_sites vector
    variable_sites = Vector{Vector{Int64}}(undef, gamma_category_count > 0 ? gamma_category_count : 1)

    # Shuffle site indices
    all_sites = shuffle(1:sequence_length)
    invariant_sites = all_sites[1:num_invariant_sites]
    variable_site_indices = all_sites[num_invariant_sites+1:end]

    if gamma_category_count == 0
        # All variable sites share the same mutation rate
        variable_sites[1] = sort(variable_site_indices)
        μ = [mutation_rate]
    else
        # Distribute variable sites into gamma categories
        base_size = div(num_variable_sites, gamma_category_count)
        remainder = rem(num_variable_sites, gamma_category_count)
        start_idx = 1

        for i in 1:gamma_category_count
            extra = i <= remainder ? 1 : 0
            end_idx = start_idx + base_size + extra - 1
            variable_sites[i] = sort(variable_site_indices[start_idx:end_idx])
            start_idx = end_idx + 1
        end

        # Compute rates using a gamma distribution
        gamma_dist = Gamma(gamma_shape, mutation_rate / (gamma_shape * num_variable_sites / sequence_length))
        quantiles = range(1/(2*gamma_category_count), step=1/gamma_category_count, length=gamma_category_count)
        μ_raw = quantile.(gamma_dist, quantiles)
        μ = μ_raw .* (mean(gamma_dist) / mean(μ_raw))
    end

    return variable_sites, μ
end