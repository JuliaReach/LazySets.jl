function convert(::Type{HPolytope}, X::LazySet)
    if !ispolytopictype(typeof(X))
        error("conversion to `HPolytope` requires a polytopic set")
    end
    return HPolytope(constraints_list(X))
end

convert(::Type{HPolytope{N,VT}}, P::HPolytope{N,VT}) where {N,VT} = P

function convert(::Type{HPolytope{N,VT}}, X::LazySet) where {N,VT}
    if !ispolytopictype(typeof(X))
        error("conversion to `HPolytope` requires a polytopic set")
    end
    return HPolytope([HalfSpace(VT(c.a), N(c.b)) for c in constraints_list(X)])
end

function load_Polyhedra_convert_HPolytope()
    return quote
        using .Polyhedra: HRep

        function convert(::Type{HPolytope}, P::HRep)
            return _convert_HPoly(HPolytope, P)
        end

        function _convert_HPoly(T, P::HRep{N}) where {N}
            VT = Polyhedra.hvectortype(P)  # NOTE: this is an internal function
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
            return T(constraints)
        end
    end
end  # load_Polyhedra_convert_HPolytope
