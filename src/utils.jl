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
    return λ, V, V⁻¹
end


"""
    mod_wrap(x::T, y::T) where T <: Integer -> T

Calculate a modified modulus that avoids returning zero. If `x` modulo `y` is zero, returns `y` instead. 
This modification ensures wrap-around behavior, making it particularly useful in scenarios such as array 
indexing or cyclic operations where a zero result is impractical.

# Arguments
- `x`: The dividend in the modulus operation.
- `y`: The divisor in the modulus operation. Must be positive to avoid division errors and logical inconsistencies.

# Returns
- The modulus of `x` by `y`, adjusted to return `y` instead of `0` if `x` is divisible by `y`.

# Example
```julia
mod_wrap(10, 5) # returns 5
mod_wrap(11, 5) # returns 1
```
# Errors
Throws an ArgumentError if y is less than or equal to zero, as modulus by non-positive numbers is not defined.
"""
function mod_wrap(x::T, y::T) where T <: Integer
    y ≤ 0 && throw(ArgumentError("The divisor 'y' must be positive. Provided y = $y"))
    x % y == 0 && return y
    return mod(x,y)
end