function rationalize(::Type{T}, X::Interval{<:AbstractFloat}, tol::Real) where {T<:Integer}
    l = rationalize(T, _min(X), tol)
    h = rationalize(T, _max(X), tol)
    return Interval(l, h)
end
