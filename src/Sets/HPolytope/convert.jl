"""
    convert(::Type{HPolytope}, X::LazySet)

Convert a polytopic set to a polytope in constraint representation.

### Input

- `HPolytope` -- target type
- `X`         -- polytopic set

### Output

The given polytope represented as a polytope in constraint representation.

### Algorithm

This method uses `constraints_list`.
"""
function convert(::Type{HPolytope}, X::LazySet)
    if !isboundedtype(typeof(X)) || !is_polyhedral(X)
        error("conversion to `HPolytope` requires a polytopic set")
    end
    return HPolytope(constraints_list(X))
end

convert(::Type{HPolytope{N,VT}}, P::HPolytope{N,VT}) where {N,VT} = P

function convert(::Type{HPolytope{N,VT}}, X::LazySet) where {N,VT}
    if !isboundedtype(typeof(X)) || !is_polyhedral(X)
        error("conversion to `HPolytope` requires a polytopic set")
    end
    return HPolytope([HalfSpace(VT(c.a), N(c.b)) for c in constraints_list(X)])
end
