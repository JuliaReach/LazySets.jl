"""
    vertices_list(cup::UnionSetArray; [apply_convex_hull]::Bool=false,
                  [backend]=nothing)

Return a list of vertices of the union of a finite number of sets.

### Input

- `cup`               -- union of a finite number of sets
- `apply_convex_hull` -- (optional, default: `false`) if `true`, post-process
                         the vertices using a convex-hull algorithm
- `backend`           -- (optional, default: `nothing`) backend for computing
                         the convex hull (see argument `apply_convex_hull`)

### Output

A list of vertices, possibly reduced to the list of vertices of the convex hull.
"""
function vertices_list(cup::UnionSetArray;
                       apply_convex_hull::Bool=false,
                       backend=nothing)
    vlist = vcat([vertices_list(Xi) for Xi in cup]...)
    if apply_convex_hull
        convex_hull!(vlist; backend=backend)
    end
    return vlist
end
