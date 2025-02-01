@commutative function distance(x::AbstractVector, L::Line; p::Real=2)
    @assert length(x) == dim(L) "incompatible dimensions $(length(x)) and $(dim(L))"

    if p != 2
        throw(ArgumentError("`distance` is only implemented for Euclidean norm"))
    end

    d = L.d  # direction of the line
    t = dot(x - L.p, d) / dot(d, d)
    return distance(x, L.p + t * d; p=p)
end
