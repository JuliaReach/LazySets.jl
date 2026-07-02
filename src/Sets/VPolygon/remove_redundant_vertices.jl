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

# see ext/LazySets/LazySetsVPolygonExt.jl
