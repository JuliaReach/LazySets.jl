function convert(::Type{HPolyhedron}, X::LazySet)
    if !ispolyhedral(X)
        throw(ArgumentError("conversion to `HPolyhedron` requires a polyhedral set"))
    end
    return HPolyhedron(constraints_list(X))
end

convert(::Type{HPolyhedron{N,VT}}, P::HPolyhedron{N,VT}) where {N,VT} = P

function convert(::Type{HPolyhedron{N,VT}}, X::LazySet) where {N,VT}
    if !ispolyhedral(X)
        throw(ArgumentError("conversion to `HPolyhedron` requires a polyhedral set"))
    end
    return HPolyhedron([HalfSpace(VT(c.a), N(c.b)) for c in constraints(X)])
end

function load_Polyhedra_convert_HPolyhedron()
    return quote
        using .Polyhedra: HRep
        using ..HPolytopeModule: _convert_HPoly

        function convert(::Type{HPolyhedron}, P::HRep)
            return _convert_HPoly(HPolyhedron, P)
        end
    end
end  # load_Polyhedra_convert_HPolyhedron
