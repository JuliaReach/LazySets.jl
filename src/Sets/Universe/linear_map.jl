"""
# Extended help

    linear_map(M::AbstractMatrix, U::Universe)

### Output

Either a `Universe` if the map is invertible, or the result of the default implementation, which
typically yields an `HPolyhedron`.
"""
@validate function linear_map(M::AbstractMatrix, U::Universe)
    if isinvertible(M)
        N = eltype(U)
        m = size(M, 1)
        return Universe{N}(m)
    end
    return _linear_map_polyhedron(M, U)
end
