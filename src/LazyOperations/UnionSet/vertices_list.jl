"""
    vertices_list(cup::UnionSet; [apply_convex_hull]::Bool=false,
                  [backend]=nothing)

Return the list of vertices of the union of two sets.

### Input

- `cup`               -- union of two sets
- `apply_convex_hull` -- (optional, default: `false`) if `true`, post-process
                         the vertices using a convex-hull algorithm
- `backend`           -- (optional, default: `nothing`) backend for computing
                         the convex hull (see argument `apply_convex_hull`)

### Output

The list of vertices, possibly reduced to the list of vertices of the convex
hull.
"""
function vertices_list(cup::UnionSet;
                       apply_convex_hull::Bool=false,
                       backend=nothing)
    vlist = vcat(vertices_list(cup.X), vertices_list(cup.Y))
    if apply_convex_hull
        convex_hull!(vlist; backend=backend)
    end
    return vlist
end
