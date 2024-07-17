"""
    σ(d::AbstractVector, P::HPolygon;
      [linear_search]::Bool=(length(P.constraints) < $BINARY_SEARCH_THRESHOLD))

Return a support vector of a polygon in a given direction.

### Input

- `d`             -- direction
- `P`             -- polygon in constraint representation
- `linear_search` -- (optional, default: see below) flag for controlling whether
                     to perform a linear search or a binary search

### Output

The support vector in the given direction.
The result is always one of the vertices; in particular, if the direction has
norm zero, any vertex is returned.

### Algorithm

Comparison of directions is performed using polar angles; see the definition of
`⪯` for two-dimensional vectors.

For polygons with $BINARY_SEARCH_THRESHOLD or more constraints we use a binary
search by default.
"""
function σ(d::AbstractVector, P::HPolygon;
           linear_search::Bool=(length(P.constraints) < BINARY_SEARCH_THRESHOLD))
    n = length(P.constraints)
    @assert n > 0 "the polygon has no constraints"

    if linear_search
        # linear search
        k = 1
        while k <= n && P.constraints[k].a ⪯ d
            k += 1
        end
    else
        # binary search
        k = binary_search_constraints(d, P.constraints)
    end

    if k == 1 || k == n + 1
        # corner cases: wrap-around in constraints list
        return element(_intersection_line2d(P.constraints[1], P.constraints[n]))
    else
        return element(_intersection_line2d(P.constraints[k], P.constraints[k - 1]))
    end
end
