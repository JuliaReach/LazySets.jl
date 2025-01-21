@commutative function distance(x::AbstractVector, L::Line; p::Real=2)
    @assert length(x) == dim(L) "a vector of length $(length(x)) is " *
                                "incompatible with a set of dimension $(dim(L))"

    d = L.d  # direction of the line
    t = dot(x - L.p, d) / dot(d, d)
    return distance(x, L.p + t * d; p=p)
end
