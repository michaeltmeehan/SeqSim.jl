
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
    κ::Float64  # Transition/Transversion ratio
    π::Vector{Float64} = fill(0.25, 4)  # Default base frequencies

    # Constructor with check for κ
    function K2P(κ::Float64, π::Vector{Float64} = fill(0.25, 4))
        κ ≤ 0 && error("The transition/transversion ratio (κ) must be positive. Received κ = $κ")
        π != fill(0.25, 4) && any(π .< 0) && error("Base frequencies must be non-negative. Received π = $π")
        π == fill(0.25, 4) && isapprox(sum(π), 1.0; atol=1e-5) || error("The sum of base frequencies (π) must be 1. Received sum = $(sum(π))")
        new(κ, π)
    end
end

K2P(κ) = K2P(κ=κ, π=π)



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
- `Q::Matrix{Float64}`: A 4x4 matrix of substitution rates, where each element [i, j] represents the rate from nucleotide i to j.
- `π::Vector{Float64}`: A vector of base frequencies for each nucleotide (A, C, G, T) that must sum to 1.

# Example
```julia
Q = [0 1 2 1; 1 0 1 2; 2 1 0 1; 1 2 1 0]
π = [0.25, 0.25, 0.25, 0.25]
model = GTR(Q, π)
"""
@with_kw mutable struct GTR{T<:Number} <: SubstitutionModel
    Q::Matrix{T}
    π::Vector{Float64}
    # Constructor with checks for rate_matrix time-reversibility
    function GTR{T}(Q::Matrix{T}, π::Vector{Float64}) where T <: Number
        sum(π) ≈ 1.0 || error("The sum of the base frequencies (π) must be 1. Received sum = $(sum(π))")
        any(π .< 0) && error("All base frequencies must be non-negative. Received frequencies = $π")
        issymmetric(Q * diagm(π)) || error("Rate matrix is not time-reversible")
    new{T}(Q, π)
    end
end


"""
    rate_matrix(π::Vector{Float64}, Q_bare::Matrix{Float64})::Matrix{Float64}

Calculate the rate matrix for a given substitution model by scaling the base rate matrix (`Q_bare`)
by the stationary distribution vector (`π`). This function multiplies each off-diagonal element
of `Q_bare` by the corresponding element in `π`, then adjusts the diagonal elements to ensure that
each row sums to zero, which is a requirement for a valid rate matrix in continuous-time Markov chains.

# Arguments
- `π::Vector{Float64}`: Stationary distribution vector, where each element represents the equilibrium frequency of a nucleotide.
- `Q_bare::Matrix{<:Number}`: A base rate matrix where each element `[i, j]` represents the unscaled rate of substitution from nucleotide `i` to `j`.

# Returns
- `Matrix{Float64}`: A valid rate matrix where each row sums to zero.

# Examples
```julia
π = [0.25, 0.25, 0.25, 0.25]  # Equal base frequencies
Q_bare = ones(4, 4) - I  # Basic rate matrix with no self-transitions
rate_matrix = rate_matrix(π, Q_bare)
"""
function rate_matrix(π::Vector{Float64}, Q_bare::Matrix{<:Number})::Matrix{Float64}
    Q = π .* Q_bare
    for i in 1:4
        Q[i,i] = -sum(Q[:,i])
    end
    return Q
end


"""
    rate_matrix(substitution_model::JC)::Matrix{Float64}

Generate the rate matrix for the Jukes-Cantor (JC) model. This model assumes equal probability
of substitution between any two different nucleotides, scaled by the stationary distribution.

# Arguments
- `substitution_model::JC`: An instance of the JC model containing base frequencies.

# Returns
- `Matrix{Float64}`: The rate matrix for the JC model, where all substitutions are equally likely.

# Examples
```julia
model = JC()
rate_matrix = rate_matrix(model)
"""
function rate_matrix(substitution_model::JC)::Matrix{Float64}
    π = substitution_model.π
    Q_bare = ones(4,4) - I
    return rate_matrix(π, Q_bare)
end


"""
    rate_matrix(substitution_model::F81)::Matrix{Float64}

Generate the rate matrix for the Felsenstein 81 (F81) model. This model allows for different
base frequencies and assumes all substitutions not involving the same base occur at rates
proportional to the equilibrium frequency of the target base.

# Arguments
- `substitution_model::F81`: An instance of the F81 model containing base frequencies.

# Returns
- `Matrix{Float64}`: The rate matrix for the F81 model, accounting for different target base frequencies.

# Examples
```julia
model = F81(π=[0.1, 0.2, 0.3, 0.4])
rate_matrix = rate_matrix(model)
"""
function rate_matrix(substitution_model::F81)::Matrix{Float64}
    π = substitution_model.π
    Q_bare = ones(4,4) - I
    return rate_matrix(π, Q_bare)
end

"""
    rate_matrix(model::K2P)::Matrix{Float64}

Calculate the rate matrix for the Kimura 2-Parameter (K2P) model, which distinguishes between transition
and transversion rates. Transitions occur at a rate κ times the rate of transversions.

# Arguments
- `model::K2P`: An instance of the K2P model containing the transition/transversion ratio (κ).

# Returns
- `Matrix{Float64}`: The rate matrix for the K2P model, where transitions are κ times more likely than transversions.

# Example
```julia
model = K2P(κ=2.0)
rate_matrix = rate_matrix(model)
"""
function rate_matrix(substitution_model::K2P)::Matrix{Float64}
    π = substitution_model.π
    κ = substitution_model.κ
    Q_bare = [0. 1. κ 1.; 1. 0. 1. κ; κ 1. 0. 1.; 1. κ 1. 0.]
    return rate_matrix(π, Q_bare)
end   


"""
    rate_matrix(substitution_model::HKY)::Matrix{Float64}

Generate the rate matrix for the Hasegawa-Kishino-Yano (HKY) model. This model differentiates
between transitions and transversions by providing a different rate for transitions (scaled by κ),
and assumes different base frequencies.

# Arguments
- `substitution_model::HKY`: An instance of the HKY model containing the transition/transversion ratio and base frequencies.

# Returns
- `Matrix{Float64}`: The rate matrix for the HKY model, differentiating transition and transversion rates.

# Examples
```julia
model = HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4])
rate_matrix = rate_matrix(model)
"""
function rate_matrix(substitution_model::HKY)::Matrix{Float64}
    π = substitution_model.π
    κ = substitution_model.κ
    Q_bare = [0. 1. κ 1.; 1. 0. 1. κ; κ 1. 0. 1.; 1. κ 1. 0.]
    return rate_matrix(π, Q_bare)
end


"""
    calculate_rate_matrix(model::GTR{T}) where T <: Number -> Matrix{T}

Calculate the rate matrix for a General Time Reversible (GTR) model. This function assumes that
the input rate matrix and the stationary probabilities already satisfy the detailed balance condition.
This function simply ensures that the rows of the resulting matrix sum to zero, which is a requirement
for a valid rate matrix in continuous-time Markov models.

# Arguments
- `model::GTR{T}`: An instance of the GTR model with predefined rate matrix and base frequencies.

# Returns
- `Matrix{T}`: A rate matrix adjusted such that each row sums to zero.

# Example
```julia
π = [0.25, 0.25, 0.25, 0.25]
rate_matrix = [0 1 2 1; 1 0 1 2; 2 1 0 1; 1 2 1 0]
model = GTR(rate_matrix, π)
adjusted_matrix = rate_matrix(model)
"""
function rate_matrix(substitution_model::GTR{<:Number})::Matrix{Float64}
    π = substitution_model.π
    Q = substitution_model.Q
    return rate_matrix(π, Q)
end