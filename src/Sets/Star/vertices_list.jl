"""
    vertices_list(X::Star; apply_convex_hull::Bool=true)

Return the list of vertices of a star.

### Input

- `X`                 -- star
- `apply_convex_hull` -- (optional, default: `true`) if `true`, apply the convex
                         hull operation to the list of vertices of the star

### Output

A list of vertices.

### Algorithm

See [`vertices_list(::LazySets.AbstractAffineMap)`](@ref).
"""
function vertices_list(X::Star; apply_convex_hull::Bool=true)
    am = convert(STAR, X)
    return vertices_list(am; apply_convex_hull=apply_convex_hull)
end
