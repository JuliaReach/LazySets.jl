"""
# Extended help

    area(V::VPolygon)

### Algorithm

See [`area(::LazySets.LazySet)`](@ref).
"""
function area(V::VPolygon)
    return _area_vlist(V.vertices; apply_convex_hull=false)
end
