
abstract type SubstitutionModel end


function rate_matrix(π::Vector{Float64}, R::Matrix{T}) where T<:Number
    Q = π .* R
    Q .-= diagm(0 => sum(Q, dims=1)[:])
    @assert all(isapprox.(sum(Q, dims=1), 0.; atol=1e-10)) "Rows of Q should sum to 0"
    @assert issymmetric(Q * Diagonal(π)) "Q should be time-reversible"
    @assert all(isapprox.(Q * π, 0.; atol=1e-10)) "Q should annihilate π"
    return Q
end


struct RateMatrixDecomposition{N}
    λ::SVector{N, Float64}
    V::SMatrix{N, N, Float64}
    V⁻¹::SMatrix{N, N, Float64}
end


function decompose(M::Matrix{<:Number})::RateMatrixDecomposition
    n = size(M, 1)
    size(M, 1) != size(M, 2) && throw(ArgumentError("Matrix must be square to perform eigen decomposition."))
    λ, V = eigen(M)
    if cond(V) > 1e12
        @warn "Matrix V is poorly conditioned (cond = $(cond(V))). Using pseudo-inverse for stability."
        V⁻¹ = pinv(V)
    else
        V⁻¹ = inv(V)
    end
    @assert isapprox(V * Diagonal(λ) * V⁻¹, M; atol=1e-10) "Eigen decomposition check failed: V * Diagonal(λ) * V⁻¹ should equal M"
    return RateMatrixDecomposition{n}(SVector{n}(λ), SMatrix{n,n}(V), SMatrix{n,n}(V⁻¹))
end


struct JC <: SubstitutionModel end

rate_matrix(substitution_model::JC) = rate_matrix(fill(0.25, 4), ones(4,4) - I)

function Base.show(io::IO, model::JC)
    print(io, "JC model (Jukes-Cantor): equal base frequencies, equal rates")
end


struct F81 <: SubstitutionModel
    π::Vector{Float64}
end

rate_matrix(substitution_model::F81) = rate_matrix(substitution_model.π, ones(4, 4) - I)

function Base.show(io::IO, model::F81)
    print(io, "F81 model: base frequencies = ", model.π)
end


struct K2P <: SubstitutionModel
    κ::Float64  # Transition/Transversion ratio
end


function rate_matrix(substitution_model::K2P)
    π = fill(0.25, 4)
    κ = substitution_model.κ
    # A/G and C/T are transitions
    R = [0 1 κ 1;
         1 0 1 κ;
         κ 1 0 1;
         1 κ 1 0]
    return rate_matrix(π, R)
end

function Base.show(io::IO, model::K2P)
    print(io, "K2P model: κ = ", model.κ, " (transition/transversion ratio)")
end


struct HKY <: SubstitutionModel
    π::Vector{Float64}
    κ::Float64
end


# Constructor with short-circuit evaluation for checking frequencies
function rate_matrix(substitution_model::HKY)
    π = substitution_model.π
    κ = substitution_model.κ
    R = [0 1 κ 1;
         1 0 1 κ;
         κ 1 0 1;
         1 κ 1 0]
    return rate_matrix(π, R)
end

function Base.show(io::IO, model::HKY)
    println(io, "HKY model:")
    println(io, "  transition/transversion ratio: ", model.κ)
    print(io, "  base frequencies: ", model.π)
end


struct GTR <: SubstitutionModel
    π::Vector{Float64}
    rates::Matrix{<:Number}
end

rate_matrix(substitution_model::GTR) = rate_matrix(substitution_model.π, substitution_model.rates)

function Base.show(io::IO, model::GTR)
    print(io, "GTR model:\n")
    print(io, "  base frequencies = ", model.π, "\n")
    print(io, "  rate matrix R =\n")
    show(io, "text/plain", model.rates)
end