"""
# Extended help

    intersection(P1::VPolygon, P2::VPolygon; apply_convex_hull::Bool=true)

### Output

A `VPolygon`, or an `EmptySet` if the intersection is empty.

### Algorithm

This function applies the [Sutherland–Hodgman polygon clipping
algorithm](https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm).
The implementation is based on the one found in
[rosetta code](http://www.rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Julia).
"""
@validate function intersection(P1::VPolygon, P2::VPolygon; apply_convex_hull::Bool=true)
    v1 = vertices_list(P1)
    v2 = vertices_list(P2)
    v12 = _intersection_vrep_2d(v1, v2)

    if isempty(v12)
        require(@__MODULE__, :LazySets; fun_name="intersection")

        N = promote_type(eltype(P1), eltype(P2))
        return EmptySet{N}(2)
    else
        return VPolygon(v12; apply_convex_hull=apply_convex_hull)
    end
end
