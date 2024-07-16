"""
    constraints_list(P::VPolygon)

Return the list of constraints defining a polygon in vertex representation.

### Input

- `P` -- polygon in vertex representation

### Output

The list of constraints of the polygon.

### Algorithm

We convert to constraint representation using `tohrep`.
"""
function constraints_list(P::VPolygon)
    return constraints_list(tohrep(P))
end
