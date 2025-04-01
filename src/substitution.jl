
abstract type SubstitutionModel end


"""
    rate_matrix(π::Vector{Float64}, R::Matrix{Float64})

Constructs a rate matrix Q given stationary distribution π and raw exchangeability matrix R.
"""
function rate_matrix(π::Vector{Float64}, R::Matrix{T}) where T<:Number
    Q = π .* R
    Q .-= diagm(0 => sum(Q, dims=1)[:])
    @assert all(isapprox.(sum(Q, dims=1), 0.; atol=1e-10)) "Rows of Q should sum to 0"
    @assert issymmetric(Q * Diagonal(π)) "Q should be time-reversible"
    @assert all(isapprox.(Q * π, 0.; atol=1e-10)) "Q should annihilate π"
    return Q
end


"""
    decompose(M::Matrix{<:Number}) -> (λ::Vector{Float64}, V::Matrix{Float64}, V⁻¹::Matrix{Float64})

Decompose a square matrix `M` of any numeric type into its eigenvalues and eigenvectors, with enhanced handling for numerical stability.

This function performs eigen decomposition on a given square matrix `M`, returning the eigenvalues as a vector (`λ`), the matrix of 
eigenvectors (`V`), and the inverse of the matrix of eigenvectors (`V⁻¹`). It handles matrices with elements of any numeric type, 
making it highly versatile. The function is robust against potential numerical instability by optionally using the pseudo-inverse 
when the eigenvector matrix is poorly conditioned.

# Arguments
- `M`: A square matrix of any numeric type. The matrix must be square to ensure a valid eigen decomposition.

# Returns
- `λ`: A vector containing the eigenvalues of `M`.
- `V`: A matrix whose columns are the eigenvectors corresponding to the eigenvalues in `λ`.
- `V⁻¹`: The inverse of the eigenvector matrix `V`, computed using the pseudo-inverse if `V` is poorly conditioned.

# Example
```julia
M = [2.0 0; 0 3.0]
λ, V, V⁻¹ = decompose(M)
```

# Notes
The function checks if M is square and will throw an ArgumentError if this condition is not met.
If the condition number of V exceeds 1e12, indicating poor conditioning, a pseudo-inverse is used instead of the regular matrix 
inverse to enhance numerical stability. This is especially useful in applications where precision and stability are crucial, such 
as in systems dynamics, stability analysis, or when dealing with ill-conditioned matrices. Diagnostic messages about the use of 
the pseudo-inverse are printed to the console for transparency and can aid in debugging and performance tuning.

# See Also
eigen: Used to compute eigenvalues and eigenvectors.
pinv: Used to compute the pseudo-inverse of a matrix.
"""
function decompose(M::Matrix{<:Number})
    size(M, 1) != size(M, 2) && throw(ArgumentError("Matrix must be square to perform eigen decomposition."))
    λ, V = eigen(M)
    if cond(V) > 1e12
        println("Matrix V is poorly condition (cond = $(cond(V))). Using pseudo-inverse for stability.")
        V⁻¹ = pinv(V)
    else
        V⁻¹ = inv(V)
    end
    @assert isapprox(V * Diagonal(λ) * V⁻¹, M; atol=1e-10) "Eigen decomposition check failed: V * Diagonal(λ) * V⁻¹ should equal M"
    return λ, V, V⁻¹
end


"""
    JC <: SubstitutionModel

The Jukes-Cantor (JC) model is a simple substitution model used in phylogenetics.
It assumes equal base frequencies and equal mutation rates between any pair of bases.
This model was originally proposed by Jukes and Cantor in their 1969 paper.

# Reference
Jukes, T. H., & Cantor, C. R. (1969). Evolution of protein molecules. In Munro, H. N. (Ed.), Mammalian Protein Metabolism (pp. 21-132). New York: Academic Press.

# Example
```julia
model = JC()
"""
# TODO: Convert diagonalized components to SMatrix{4, 4, Float64} from StaticArrays.jl for performance
struct JC <: SubstitutionModel
    π::Vector{Float64}
    Q::Matrix{Float64}
    λ::SVector{4, Float64}
    V::SMatrix{4, 4, Float64}
    V⁻¹::SMatrix{4, 4, Float64}
end


function JC()
    π = fill(0.25, 4)
    R = ones(4, 4) - I
    Q = rate_matrix(π, R)
    λ, V, V⁻¹ = decompose(Q)
    return JC(π, Q, SVector{4}(λ), SMatrix{4,4}(V), SMatrix{4,4}(V⁻¹))
end


"""
    F81 <: SubstitutionModel

The Felsenstein 81 (F81) model extends the Jukes-Cantor model by allowing different
base frequencies but maintaining the assumption of equal mutation rates across different pairs.
This model was introduced by Joseph Felsenstein in 1981, providing a method to consider
different nucleotide frequencies in evolutionary studies.

# Reference
Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. 
Journal of Molecular Evolution, 17(6), 368-376.

# Example
```julia
model = F81([0.1, 0.2, 0.3, 0.4])
"""
struct F81 <: SubstitutionModel
    π::Vector{Float64}
    Q::Matrix{Float64}
    λ::SVector{4, Float64}
    V::SMatrix{4, 4, Float64}
    V⁻¹::SMatrix{4, 4, Float64}
end


function F81(π::Vector{Float64})
    any(π .< 0) && throw(ArgumentError("All frequencies must be non-negative. Received frequencies = $π"))
    sum(π) ≈ 1.0 || throw(ArgumentError("Base frequencies must sum to 1"))
    R = ones(4, 4) - I
    Q = rate_matrix(π, R)
    λ, V, V⁻¹ = decompose(Q)
    return F81(π, Q, SVector{4}(λ), SMatrix{4,4}(V), SMatrix{4,4}(V⁻¹))
end


"""
    K2P <: SubstitutionModel

The Kimura 2-Parameter (K2P) model, also known as the Kimura two-parameter model, is a substitution model used in phylogenetics 
for nucleotide sequences. It distinguishes between transition and transversion mutation rates, assuming two types of 
substitution rates but equal base frequencies. This model provides a more realistic depiction of molecular evolution than 
the simpler Jukes-Cantor model by accounting for the fact that transitions often occur at different rates than transversions.

# Reference
Kimura, M. (1980). A simple method for estimating evolutionary rates of base substitutions through comparative 
studies of nucleotide sequences. Journal of Molecular Evolution, 16(2), 111-120.

# Example
```julia
model = K2P(2.0)
"""
struct K2P <: SubstitutionModel
    π::Vector{Float64}
    κ::Float64  # Transition/Transversion ratio
    Q::Matrix{Float64}
    λ::SVector{4, Float64}
    V::SMatrix{4, 4, Float64}
    V⁻¹::SMatrix{4, 4, Float64}
end


function K2P(κ::Float64)
    κ > 0 || throw(ArgumentError("κ must be positive"))
    π = fill(0.25, 4)
    # A/G and C/T are transitions
    R = [0 1 κ 1;
         1 0 1 κ;
         κ 1 0 1;
         1 κ 1 0]
    Q = rate_matrix(π, R)
    λ, V, V⁻¹ = decompose(Q)
    return K2P(π, κ, Q, SVector{4}(λ), SMatrix{4,4}(V), SMatrix{4,4}(V⁻¹))
end


"""
    HKY <: SubstitutionModel

The Hasegawa-Kishino-Yano (HKY) model introduces a transition/transversion ratio,
providing a more realistic representation of molecular evolution by allowing different rates for transitions and transversions.
This model was proposed by Hasegawa, Kishino, and Yano in 1985 and is particularly useful
for models involving nucleotide sequences where transitions occur more frequently than transversions.

# Reference
Hasegawa, M., Kishino, H., & Yano, T. (1985). Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. 
Journal of Molecular Evolution, 22(2), 160-174.

# Example
```julia
model = HKY([0.1, 0.2, 0.3, 0.4], 2.0)
"""
struct HKY <: SubstitutionModel
    π::Vector{Float64}
    κ::Float64
    Q::Matrix{Float64}
    λ::SVector{4, Float64}
    V::SMatrix{4, 4, Float64}
    V⁻¹::SMatrix{4, 4, Float64}
end


# Constructor with short-circuit evaluation for checking frequencies
function HKY(π::Vector{Float64}, κ::Float64)
    any(π .< 0) && throw(ArgumentError("All frequencies must be non-negative. Received frequencies = $π"))
    sum(π) ≈ 1.0 || throw(ArgumentError("The sum of the frequencies must be 1. Received sum = $(sum(π))"))
    κ ≤ 0 && throw(ArgumentError("The transition/transversion ratio (κ) must be positive. Received κ = $κ"))
    R = [0 1 κ 1;
         1 0 1 κ;
         κ 1 0 1;
         1 κ 1 0]
    Q = rate_matrix(π, R)
    λ, V, V⁻¹ = decompose(Q)
    return HKY(π, κ, Q, SVector{4}(λ), SMatrix{4,4}(V), SMatrix{4,4}(V⁻¹))
end


"""
    GTR <: SubstitutionModel

The General Time Reversible (GTR) model is the most comprehensive model for nucleotide substitution,
allowing different rates for all types of substitutions and different base frequencies.
This model is essential for accurately modeling the complex evolutionary patterns in DNA sequences.

# Example
```julia
rates = [0 1 2 1; 1 0 1 2; 2 1 0 1; 1 2 1 0]
π = [0.25, 0.25, 0.25, 0.25]
model = GTR(π, rates)
"""
struct GTR <: SubstitutionModel
    π::Vector{Float64}
    rates::Matrix{<:Number}
    Q::Matrix{Float64}
    λ::SVector{4, Float64}
    V::SMatrix{4, 4, Float64}
    V⁻¹::SMatrix{4, 4, Float64}
end


# Constructor with checks for rate_matrix time-reversibility
function GTR(π::Vector{Float64}, rates::Matrix{T}) where T<:Number
    length(π) != 4 && throw(ArgumentError("Base frequencies (π) must be a vector of length 4. Received length = $(length(π))"))
    sum(π) ≈ 1.0 || throw(ArgumentError("The sum of the base frequencies (π) must be 1. Received sum = $(sum(π))"))
    any(π .< 0) && throw(ArgumentError("All base frequencies must be non-negative. Received frequencies = $π"))
    size(rates) != (4, 4) && throw(ArgumentError("Rate matrix must be 4x4. Received size = $(size(rates))"))
    issymmetric(rates * diagm(π)) || throw(ArgumentError("Rate matrix is not time-reversible"))
    Q = rate_matrix(π, rates)
    λ, V, V⁻¹ = decompose(Q)
    return GTR(π, rates, Q, SVector{4}(λ), SMatrix{4,4}(V), SMatrix{4,4}(V⁻¹))
end


function Base.show(io::IO, model::SubstitutionModel)
    model_name = nameof(typeof(model))
    println(io, "$model_name Substitution Model")

    println(io, "π = ", model.π)

    if model isa HKY || model isa K2P
        println(io, "κ = ", model.κ)
    end

    if model isa GTR
        println(io, "Rates = ")
        show(io, "text/plain", model.rates)
    end
end
