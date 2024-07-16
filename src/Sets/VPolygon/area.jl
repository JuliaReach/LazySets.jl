"""
    area(V::VPolygon)

Compute the area of a polygon in vertex representation.

### Input

- `V` -- polygon in vertex representation

### Output

A number representing the area of `V`.

### Algorithm

See [`area(::LazySets.LazySet)`](@ref).
"""
function area(V::VPolygon)
    return _area_vlist(V.vertices; apply_convex_hull=false)
end
