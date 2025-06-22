function rationalize(::Type{T}, X::Interval{<:AbstractFloat}, tol::Real) where {T<:Integer}
    l = rationalize(T, min(X), tol)
    h = rationalize(T, max(X), tol)
    return Interval(l, h)
end
