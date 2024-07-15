"""
    remove_redundant_vertices(P::VPolygon;
                              [algorithm]::String="monotone_chain")

Return a polygon obtained by removing the redundant vertices of the given
polygon.

### Input

- `P`         -- polygon in vertex representation
- `algorithm` -- (optional, default: "monotone_chain") the algorithm used to
                 compute the convex hull

### Output

A new polygon such that its vertices are the convex hull of the given polygon.

### Algorithm

See [`remove_redundant_vertices!(::VPolygon)`](@ref).
"""
function remove_redundant_vertices(P::VPolygon;
                                   algorithm::String="monotone_chain")
    return remove_redundant_vertices!(copy(P); algorithm=algorithm)
end

"""
    remove_redundant_vertices!(P::VPolygon;
                               [algorithm]::String="monotone_chain")

Remove the redundant vertices from the given polygon in-place.

### Input

- `P`         -- polygon in vertex representation
- `algorithm` -- (optional, default: "monotone_chain") the algorithm used to
                 compute the convex hull

### Output

The modified polygon whose redundant vertices have been removed.

### Algorithm

A convex-hull algorithm is used to compute the convex hull of the vertices of
the polygon `P`; see `?convex_hull` for details on the available algorithms.
The vertices are sorted in counter-clockwise fashion.
"""
function remove_redundant_vertices!(P::VPolygon;
                                    algorithm::String="monotone_chain")
    convex_hull!(P.vertices; algorithm=algorithm)
    return P
end
