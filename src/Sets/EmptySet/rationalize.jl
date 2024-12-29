function rationalize(::Type{T}, ∅::EmptySet{<:AbstractFloat}, tol::Real) where {T<:Integer}
    return EmptySet{Rational{T}}(dim(∅))
end
