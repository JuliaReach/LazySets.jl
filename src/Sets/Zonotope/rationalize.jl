function rationalize(::Type{T}, Z::Zonotope{<:AbstractFloat}, tol::Real) where {T<:Integer}
    c = rationalize(T, Z.center, tol)
    G = rationalize(T, Z.generators, tol)
    return Zonotope(c, G)
end
