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

function load_Polyhedra_convert_HPolytope()
    return quote
        using .Polyhedra: HRep

        function convert(::Type{HPolytope}, P::HRep{N}) where {N}
            VT = Polyhedra.hvectortype(P)
            constraints = Vector{HalfSpace{N,VT}}()
            for hi in Polyhedra.allhalfspaces(P)
                a, b = hi.a, hi.β
                if isapproxzero(norm(a))
                    @assert b >= zero(N) "the half-space is inconsistent since it has a zero " *
                                         "normal direction but the constraint is negative"
                    continue
                end
                push!(constraints, HalfSpace(a, b))
            end
            return HPolytope(constraints)
        end
    end
end  # load_Polyhedra_convert_HPolytope
