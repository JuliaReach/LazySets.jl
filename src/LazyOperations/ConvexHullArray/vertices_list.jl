"""
    vertices_list(cha::ConvexHullArray; [apply_convex_hull]::Bool=true,
                  [backend]=nothing, [prune]::Bool=apply_convex_hull)

Return a list of vertices of the convex hull of a finite number of sets.

### Input

- `cha`               -- convex hull of a finite number of sets
- `apply_convex_hull` -- (optional, default: `true`) if `true`, post-process the
                         vertices using a convex-hull algorithm
- `backend`           -- (optional, default: `nothing`) backend for computing
                         the convex hull (see argument `apply_convex_hull`)
- `prune`             -- (optional, default: `apply_convex_hull`) alias for
                         `apply_convex_hull`

### Output

A list of vertices.
"""
@validate function vertices_list(cha::ConvexHullArray;
                                 apply_convex_hull::Bool=true,
                                 backend=nothing,
                                 prune::Bool=apply_convex_hull)
    vlist = vcat([vertices_list(Xi) for Xi in cha]...)
    if apply_convex_hull || prune
        convex_hull!(vlist; backend=backend)
    end
    return vlist
end
