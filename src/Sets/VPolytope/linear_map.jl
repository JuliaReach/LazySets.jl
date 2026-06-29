"""
# Extended help

    linear_map(M::AbstractMatrix, P::VPolytope; [apply_convex_hull]::Bool=false)

### Input

- `apply_convex_hull` -- (optional, default: `false`) flag for applying a convex
                         hull to eliminate redundant vertices

### Algorithm

The linear map ``M`` is applied to each vertex of the given set ``P``, obtaining
a polytope in vertex representation. The output type is again a `VPolytope`.
"""
@validate function linear_map(M::AbstractMatrix, P::VPolytope; apply_convex_hull::Bool=false)
    return _linear_map_vrep(M, P; apply_convex_hull=apply_convex_hull)
end

# see ext/LazySets/LazySetsVPolytopeExt.jl
function _linear_map_vrep(M::AbstractMatrix, P,
                          algo::LinearMapVRep=LinearMapVRep(nothing);  # ignored
                          apply_convex_hull::Bool=false)
    return error()
end
