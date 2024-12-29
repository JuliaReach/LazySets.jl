@commutative function distance(x::AbstractVector, L::Line; p::Real=2)
    d = L.d  # direction of the line
    t = dot(x - L.p, d) / dot(d, d)
    return distance(x, L.p + t * d; p=p)
end
