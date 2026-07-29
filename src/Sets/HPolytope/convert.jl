function convert(::Type{HPolytope}, X::LazySet)
    @assert ispolytopictype(typeof(X)) "conversion to `HPolytope` requires a polytopic set"

    return HPolytope(constraints_list(X))
end

convert(::Type{HPolytope{N,VT}}, P::HPolytope{N,VT}) where {N,VT} = P

function convert(::Type{HPolytope{N,VT}}, X::LazySet) where {N,VT}
    @assert ispolytopictype(typeof(X)) "conversion to `HPolytope` requires a polytopic set"

    return HPolytope([HalfSpace(VT(c.a), N(c.b)) for c in constraints_list(X)])
end
