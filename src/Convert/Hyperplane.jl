function Base.convert(::Type{Hyperplane}, P::HPolyhedron; skip_check::Bool=false)
    # check that the number of constraints is fine
    if !skip_check && !ishyperplanar(P)
        throw(ArgumentError("the polyhedron is not hyperplanar: $P"))
    end

    # construct hyperplane from first constraint
    c1 = @inbounds first(constraints_list(P))
    return Hyperplane(c1.a, c1.b)
end

function Base.convert(::Type{Hyperplane}, L::Line)
    if dim(L) != 2
        throw(ArgumentError("can only convert a 2D Line to a Hyperplane"))
    end

    a = @inbounds [L.d[2], -L.d[1]]  # rotate direction by 90°
    b = dot(a, L.p)  # substitute line's point
    return Hyperplane(a, b)
end
