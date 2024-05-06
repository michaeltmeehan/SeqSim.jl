mutable struct SiteModel
    mutation_rate::Float64
    gamma_category_count::Int
    gamma_shape::Float64
    proportion_invariant::Float64
    substitution_model::SubstitutionModel
end