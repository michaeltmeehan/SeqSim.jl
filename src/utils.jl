
function decompose(M::Matrix{Float64})
    U, V = eigen(M)
    D = diagm(U)
    V⁻¹ = inv(V)
    return D, V, V⁻¹
end


function mod2(x::T, y::T) where T <: Integer
    x % y == 0 && return y
    return mod(x,y)
end