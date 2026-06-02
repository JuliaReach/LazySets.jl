"""
    linear_map(M::AbstractMatrix, cpa::CartesianProductArray)

Concrete linear map of a Cartesian product of a finite number of (polyhedral)
sets.

### Input

- `M`   -- matrix
- `cpa` -- Cartesian product of a finite number of sets

### Output

A polyhedron or polytope.
"""
@validate function linear_map(M::AbstractMatrix, cpa::CartesianProductArray)
    return _linear_map_cartesian_product(M, cpa)
end
