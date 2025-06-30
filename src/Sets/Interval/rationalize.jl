function rationalize(::Type{T}, X::Interval{<:AbstractFloat}, tol::Real) where {T<:Integer}
    l = rationalize(T, X.dat.lo, tol)
    h = rationalize(T, X.dat.hi, tol)
    return Interval(l, h)
end
