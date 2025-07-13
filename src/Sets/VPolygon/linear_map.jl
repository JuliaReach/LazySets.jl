"""
    linear_map(M::AbstractMatrix, P::VPolygon; [apply_convex_hull]::Bool=false)

Concrete linear map of a polygon in vertex representation.

### Input

- `M`                 -- matrix
- `P`                 -- polygon in vertex representation
- `apply_convex_hull` -- (optional; default: `false`) flag to apply a
                         convex-hull operation (only relevant for
                         higher-dimensional maps)

### Output

The type of the result depends on the dimension. in 1D it is an interval, in 2D
it is a `VPolygon`, and in all other cases it is a `VPolytope`.

### Algorithm

This implementation uses the internal `_linear_map_vrep` method.
"""
@validate function linear_map(M::AbstractMatrix, P::VPolygon; apply_convex_hull::Bool=false)
    return _linear_map_vrep(M, P; apply_convex_hull=apply_convex_hull)
end
