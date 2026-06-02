"""
    vertices_list(ch::ConvexHull; [apply_convex_hull]::Bool=true,
                  [backend]=nothing)

Return a list of vertices of the convex hull of two sets.

### Input

- `ch`                -- convex hull of two sets
- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         vertices using a convex-hull algorithm
- `backend`           -- (optional, default: `nothing`) backend for computing
                         the convex hull (see argument `apply_convex_hull`)

### Output

A list of vertices.
"""
@validate function vertices_list(ch::ConvexHull;
                                 apply_convex_hull::Bool=true,
                                 backend=nothing)
    vlist = vcat(vertices_list(ch.X), vertices_list(ch.Y))
    if apply_convex_hull
        convex_hull!(vlist; backend=backend)
    end
    return vlist
end
