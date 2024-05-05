
abstract type SubstitutionModel end


"""
    JC <: SubstitutionModel

The Jukes-Cantor (JC) model is a simple substitution model used in phylogenetics.
It assumes equal base frequencies and equal mutation rates between any pair of bases.
This model was originally proposed by Jukes and Cantor in their 1969 paper.

# Fields
- `π::Vector{Float64}`: A vector of base frequencies, defaulting to `[0.25, 0.25, 0.25, 0.25]`, representing equal probabilities for each nucleotide (A, C, G, T).

# Reference
Jukes, T. H., & Cantor, C. R. (1969). Evolution of protein molecules. In Munro, H. N. (Ed.), Mammalian Protein Metabolism (pp. 21-132). New York: Academic Press.

# Example
```julia
model = JC()
"""
@with_kw mutable struct JC <: SubstitutionModel
    π::Vector{Float64}=fill(0.25, 4)
end


"""
    F81 <: SubstitutionModel

The Felsenstein 81 (F81) model extends the Jukes-Cantor model by allowing different
base frequencies but maintaining the assumption of equal mutation rates across different pairs.
This model was introduced by Joseph Felsenstein in 1981, providing a method to consider
different nucleotide frequencies in evolutionary studies.

# Fields
- `π::Vector{Float64}`: A vector of base frequencies that should sum to 1.

# Reference
Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. 
Journal of Molecular Evolution, 17(6), 368-376.

# Example
```julia
model = F81(π=[0.1, 0.2, 0.3, 0.4])
"""
@with_kw mutable struct F81 <: SubstitutionModel
    π::Vector{Float64}

    # Adding a constructor to check the sum and positivity of frequencies
    function F81(π::Vector{Float64})
        any(π .< 0) && error("All frequencies must be non-negative. Received frequencies = $π")
        sum(π) ≈ 1.0 || error("The sum of the frequencies must be 1. Received sum = $(sum(π))")
        new(π)
    end
end


"""
    K2P <: SubstitutionModel

The Kimura 2-Parameter (K2P) model, also known as the Kimura two-parameter model, is a substitution model used in phylogenetics 
for nucleotide sequences. It distinguishes between transition and transversion mutation rates, assuming two types of 
substitution rates but equal base frequencies. This model provides a more realistic depiction of molecular evolution than 
the simpler Jukes-Cantor model by accounting for the fact that transitions often occur at different rates than transversions.

# Fields
- `κ::Float64`: The transition/transversion ratio, a measure of the relative likelihood of transitions compared to transversions.

# Reference
Kimura, M. (1980). A simple method for estimating evolutionary rates of base substitutions through comparative 
studies of nucleotide sequences. Journal of Molecular Evolution, 16(2), 111-120.

# Example
```julia
model = K2P(κ=2.0)
"""
@with_kw mutable struct K2P <: SubstitutionModel
    κ::Float64 # Transition/Transversion ratio
    # Constructor with check for κ
    function K2P(κ::Float64)
        κ ≤ 0 && error("The transition/transversion ratio (κ) must be positive. Received κ = $κ")
        new(κ)
    end
end


"""
    HKY <: SubstitutionModel

The Hasegawa-Kishino-Yano (HKY) model introduces a transition/transversion ratio,
providing a more realistic representation of molecular evolution by allowing different rates for transitions and transversions.
This model was proposed by Hasegawa, Kishino, and Yano in 1985 and is particularly useful
for models involving nucleotide sequences where transitions occur more frequently than transversions.

# Fields
- `κ::Float64`: The transition/transversion ratio, a measure of the relative likelihood of transition to transversion mutations.
- `π::Vector{Float64}`: A vector of base frequencies that should sum to 1.

# Reference
Hasegawa, M., Kishino, H., & Yano, T. (1985). Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. 
Journal of Molecular Evolution, 22(2), 160-174.

# Example
```julia
model = HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4])
"""
@with_kw mutable struct HKY <: SubstitutionModel
    κ::Float64
    π::Vector{Float64}

    # Constructor with short-circuit evaluation for checking frequencies
    function HKY(κ::Float64, π::Vector{Float64})
        any(π .< 0) && error("All frequencies must be non-negative. Received frequencies = $π")
        sum(π) ≈ 1.0 || error("The sum of the frequencies must be 1. Received sum = $(sum(π))")
        κ ≤ 0 && error("The transition/transversion ratio (κ) must be positive. Received κ = $κ")
        new(κ, π)
    end
end


"""
    GTR <: SubstitutionModel

The General Time Reversible (GTR) model is the most comprehensive model for nucleotide substitution,
allowing different rates for all types of substitutions and different base frequencies.
This model is essential for accurately modeling the complex evolutionary patterns in DNA sequences.

# Fields
- `rate_matrix::Matrix{Float64}`: A 4x4 matrix of substitution rates, where each element [i, j] represents the rate from nucleotide i to j.
- `π::Vector{Float64}`: A vector of base frequencies for each nucleotide (A, C, G, T) that must sum to 1.

# Example
```julia
rate_matrix = [0 1 2 1; 1 0 1 2; 2 1 0 1; 1 2 1 0]
π = [0.25, 0.25, 0.25, 0.25]
model = GTR(rate_matrix, π)
"""
@with_kw mutable struct GTR <: SubstitutionModel
    rate_matrix::Matrix{Float64}
    π::Vector{Float64}
    # Constructor with checks for rate_matrix time-reversibility
    function GTR(rate_matrix::Matrix{Float64}, π::Vector{Float64})
        sum(π) ≈ 1.0 || error("The sum of the base frequencies (π) must be 1. Received sum = $(sum(π))")
        any(π .< 0) && error("All base frequencies must be non-negative. Received frequencies = $π")
    
        # Check for time-reversibility
        for i in 1:length(π)
            for j in 1:length(π)
                if i != j && abs(π[i] * rate_matrix[i, j] - π[j] * rate_matrix[j, i]) > 1e-5
                    error("Rate matrix is not time-reversible. Failed at (i=$i, j=$j)")
                end
            end
        end
    new(rate_matrix, π)
    end
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
    κ = substitution_model.κ
    Q_bare = [0. 1. κ 1.; 1. 0. 1. κ; κ 1. 0. 1.; 1. κ 1. 0.]
    return rate_matrix(π, Q_bare)
end