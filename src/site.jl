"""
    SiteModel

A `SiteModel` encapsulates the evolutionary model parameters for each site in a 
phylogenetic analysis. This structure is designed to handle variability across 
sites with respect to mutation rates, rate heterogeneity (gamma-distributed), 
proportion of invariant sites, and the type of nucleotide substitution model 
employed.

# Fields
- `mutation_rate::Float64`: The overall mutation rate at each site. Must be positive.
- `gamma_category_count::Int`: The number of rate categories for the gamma distribution describing rate heterogeneity across sites. Must be at least 1.
- `gamma_shape::Float64`: The shape parameter of the gamma distribution for rate heterogeneity. Must be positive.
- `proportion_invariant::Float64`: The proportion of sites that are invariant (no change across the phylogeny). Must be between 0 and 1.
- `substitution_model::SubstitutionModel`: The substitution model to be used. This can be any model derived from the `SubstitutionModel` abstract base class, allowing for flexible adaptation to different evolutionary scenarios.

# Constructor
The constructor checks for the validity of all input parameters to prevent erroneous or unphysical values from being used in the model.

# Example
```julia
model = SiteModel(
    mutation_rate = 1.0,
    gamma_category_count = 4,
    gamma_shape = 0.5,
    proportion_invariant = 0.1,
    substitution_model = JC()
)
"""
@with_kw mutable struct SiteModel
    mutation_rate::Float64 # Overall mutation rate at each site
    gamma_category_count::Int # Number of gamma rate categories
    gamma_shape::Float64 # Shape parameter for the gamma distribution
    proportion_invariant::Float64 # Proportion of invariant sites
    substitution_model::SubstitutionModel # Substitution model used

    # Constructor with validation
    function SiteModel(
        mutation_rate::Float64,
        gamma_category_count::Int,
        gamma_shape::Float64,
        proportion_invariant::Float64,
        substitution_model::SubstitutionModel
        )
        mutation_rate <= 0. && error("Mutation rate must be positive.")
        gamma_category_count < 1 && error("Gamma category count must be at least 1.")
        gamma_shape <= 0. && error("Gamma shape must be positive.")
        (proportion_invariant < 0. || proportion_invariant > 1.) && error("Proportion invariant must be between 0 and 1.")
    
    new(mutation_rate, gamma_category_count, gamma_shape, proportion_invariant, substitution_model)
    end
end