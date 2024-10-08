"""
    convert(::Type{HPolyhedron}, X::LazySet)

Convert a polyhedral set to a polyhedron in constraint representation.

### Input

- `HPolyhedron` -- target type
- `X`           -- polyhedral set

### Output

The given set represented as a polyhedron in constraint representation.

### Algorithm

This method uses `constraints_list`.
"""
function convert(::Type{HPolyhedron}, X::LazySet)
    if !ispolyhedral(X)
        error("conversion to `HPolyhedron` requires a polyhedral set")
    end
    return HPolyhedron(constraints_list(X))
end

convert(::Type{HPolyhedron{N,VT}}, P::HPolyhedron{N,VT}) where {N,VT} = P

function convert(::Type{HPolyhedron{N,VT}}, X::LazySet) where {N,VT}
    if !ispolyhedral(X)
        error("conversion to `HPolyhedron` requires a polyhedral set")
    end
    return HPolyhedron([HalfSpace(VT(c.a), N(c.b)) for c in constraints(X)])
end

function load_Polyhedra_convert_HPolyhedron()
    return quote
        using .Polyhedra: HRep
        using ..HPolytopeModule: _convert_HPoly

        """
             convert(::Type{HPolyhedron}, P::HRep)

        Convert an `HRep` polyhedron from `Polyhedra.jl` to a polyhedron in constraint
        representation .

        ### Input

        - `HPolyhedron` -- target type
        - `P`           -- `HRep` polyhedron

        ### Output

        An `HPolyhedron`.
        """
        function convert(::Type{HPolyhedron}, P::HRep)
            return _convert_HPoly(HPolyhedron, P)
        end
    end
end  # load_Polyhedra_convert_HPolyhedron
