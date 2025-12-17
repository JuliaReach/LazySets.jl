"""
# Extended help

    vertices_list(X::Star; apply_convex_hull::Bool=true)

### Algorithm

See [`vertices_list(::LazySets.AbstractAffineMap)`](@ref).
"""
@validate function vertices_list(X::Star; apply_convex_hull::Bool=true)
    am = convert(STAR, X)
    return vertices_list(am; apply_convex_hull=apply_convex_hull)
end
