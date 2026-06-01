"""
    linear_map(M::AbstractMatrix, cp::CartesianProduct)

Concrete linear map of a (polyhedral) Cartesian product.

### Input

- `M`  -- matrix
- `cp` -- Cartesian product

### Output

A polytope if `cp` is bounded and a polyhedron otherwise.

### Algorithm

We convert the Cartesian product to constraint representation and then call
`linear_map` on the corresponding polyhedron.

This is a fallback implementation and will fail if the wrapped sets are not
polyhedral.
"""
@validate function linear_map(M::AbstractMatrix, cp::CartesianProduct)
    return _linear_map_cartesian_product(M, cp)
end

function _linear_map_cartesian_product(M, cp)
    # use constraint representation
    T = isbounded(cp) ? HPolytope : HPolyhedron
    P = T(constraints_list(cp))
    return linear_map(M, P)
end
