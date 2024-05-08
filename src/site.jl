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
    gamma_category_count::Int64=0   # Number of bins for discretized site-specific mutation rate distribution
    gamma_shape::Float64=0. # Shape parameter for the gamma distribution
    proportion_invariant::Float64=0. # Proportion of invariant sites
    substitution_model::SubstitutionModel=JC() # Substitution model used

    # Constructor with validation
    function SiteModel(
        mutation_rate::Float64,
        gamma_category_count::Int64=0,
        gamma_shape::Float64=0.,
        proportion_invariant::Float64=0.,
        substitution_model::SubstitutionModel=JC()
        )
        mutation_rate <= 0. && error("Mutation rate must be positive.")
        gamma_shape < 0 && error("Gamma category count must be non-negative.")
        gamma_shape < 0. && error("Gamma shape must be non-negative.")
        (proportion_invariant < 0. || proportion_invariant ≥ 1.) && error("Proportion invariant must be between 0 and 1.")
        
    new(mutation_rate, gamma_category_count, gamma_shape, proportion_invariant, substitution_model)
    end
end


"""
    site_rates(model::SiteModel, n::Int64)::Vector{Float64}

Calculate site-specific mutation rates for a given number of sites (`n`) based on a specified `SiteModel`.

The function generates mutation rates using a gamma distribution for variable sites, adjusted according to 
the proportion of invariant sites specified in the model. If `gamma_shape` is set to zero, all sites are 
assigned the `mutation_rate` directly.

# Arguments
- `model::SiteModel`: The site model containing parameters for mutation rate calculations. It should include:
  - `mutation_rate`: The average mutation rate when sites are variable.
  - `gamma_shape`: The shape parameter of the gamma distribution used for variable sites.
  - `proportion_invariant`: The proportion of sites that are invariant (no mutation).
  - `substitution_model`: The substitution model (currently not used in the rate calculations).
- `n::Int64`: The number of sites for which to generate mutation rates.

# Returns
- `Vector{Float64}`: A vector of mutation rates for each site. Sites designated as invariant will have a mutation rate of 0.

# Example
```julia
model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, proportion_invariant=0.1, substitution_model=JC())
rates = site_rates(model, 100)  # Generate rates for 100 sites
```
"""
function site_rates(model::SiteModel, n::Int64)::Vector{Float64}

    n ≤ 0 && throw(ArgumentError("n must be a positive integer. Provided n = $n"))

    @unpack mutation_rate, gamma_category_count, gamma_shape, proportion_invariant, substitution_model = model
    rates = fill(0., n)
    d = Gamma(gamma_shape, mutation_rate / (gamma_shape * (1. - proportion_invariant)))
    q = quantile(d, range(start=1. /(2. * gamma_category_count), step=1. / gamma_category_count, length = gamma_category_count))
    variable_sites = rand(n) .> proportion_invariant
    rates[variable_sites] = [rand(d_disc) for _ in 1:count(variable_sites)]
    return rates
end


function rate_cat(model::SiteModel, n::Int64)::Vector{Int64}
    n ≤ 0 && throw(ArgumentError("n must be a positive integer. Provided n = $n"))
    @unpack mutation_rate, gamma_category_count, gamma_shape, proportion_invariant, substitution_model = model
    rates = fill(1, n)
    variable_sites = rand(n) .> proportion_invariant
    rates[variable_sites] = rand(2:(1 + gamma_category_count), count(variable_sites))
    return rates
end



function discretize_dist(d::T, n::Int64) where T <: Distribution{Univariate, Continuous}
    n ≤ 0 && throw(ArgumentError("n must be a positive integer. Provided n = $n"))
    locations = quantile(d, range(start=1. /(2. *n), step=1. / n, length = n))
    return MixtureModel([Dirac(loc) for loc in locations])
end