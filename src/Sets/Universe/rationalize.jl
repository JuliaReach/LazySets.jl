function rationalize(::Type{T}, U::Universe{<:AbstractFloat}, tol::Real) where {T<:Integer}
    return Universe{Rational{T}}(dim(U))
end
