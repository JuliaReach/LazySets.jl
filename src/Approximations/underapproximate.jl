"""
    underapproximate(X::LazySet{N}, dirs::AbstractDirections;
                    [apply_convex_hull]::Bool=false) where {N<:Real}

Compute the underapproximation of a convex set by sampling support vectors.

### Input

- `X`                 -- set
- `dirs`              -- directions
- `apply_convex_hull` -- (optional, default: `false`) if `true`, post-process
                         the support vectors with a convex hull operation
### Output

The `VPolytope` obtained by taking the convex hull of the support vectors of `X`
along the directions determined by `dirs`.
"""
function underapproximate(X::LazySet{N}, dirs::AbstractDirections;
                          apply_convex_hull::Bool=false) where {N<:Real}
    @assert LazySets.dim(X) == LazySets.dim(dirs)
    vinner = Vector{Vector{N}}(undef, length(dirs)) # TODO could be generalized to eltype(dirs)
    @inbounds for (i, di) in enumerate(dirs)
        vinner[i] = σ(Vector(di), X)
    end
    if apply_convex_hull
        convex_hull!(vinner)
    end
    return VPolytope(vinner)
end
