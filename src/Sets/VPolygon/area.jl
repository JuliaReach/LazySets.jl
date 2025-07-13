"""
# Extended help

    area(V::VPolygon)

### Algorithm

See [`area(::LazySets.LazySet)`](@ref).
"""
@validate function area(V::VPolygon)
    return _area_vlist_2D(V.vertices; apply_convex_hull=false)
end
