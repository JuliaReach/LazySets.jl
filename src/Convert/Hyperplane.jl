function Base.convert(::Type{Hyperplane}, P::HPolyhedron; skip_check::Bool=false)
    # check that the number of constraints is fine
    if !skip_check && !ishyperplanar(P)
        throw(ArgumentError("the polyhedron is not hyperplanar: $P"))
    end

    # construct hyperplane from first constraint
    c1 = @inbounds first(constraints_list(P))
    return Hyperplane(c1.a, c1.b)
end
