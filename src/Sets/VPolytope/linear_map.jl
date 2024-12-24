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
function linear_map(M::AbstractMatrix, P::VPolytope;
                    apply_convex_hull::Bool=false)
    @assert dim(P) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(P))"

    return _linear_map_vrep(M, P; apply_convex_hull=apply_convex_hull)
end

@inline function _linear_map_vrep(M::AbstractMatrix, P::VPolytope,
                                  algo::LinearMapVRep=LinearMapVRep(nothing);  # ignored
                                  apply_convex_hull::Bool=false)
    vlist = broadcast(v -> M * v, P.vertices)
    if apply_convex_hull
        convex_hull!(vlist)
    end
    return VPolytope(vlist)
end
