"""
    underapproximate(X::LazySet{N}, dirs::AbstractDirections;
                    [apply_convex_hull]::Bool=false) where {N}

Compute the underapproximation of a convex set by sampling support vectors.

### Input

- `X`                 -- set
- `dirs`              -- directions
- `apply_convex_hull` -- (optional, default: `false`) if `true`, post-process
                         the support vectors with a convex hull operation
### Output

The `VPolytope` obtained by taking the convex hull of the support vectors of `X`
along the directions determined by `dirs`.

### Notes

Since the support vectors are not always unique, this algorithm may return
a strict underapproximation even if the set can be exactly approximated using
the given template.
"""
function underapproximate(X::LazySet{N}, dirs::AbstractDirections;
                          apply_convex_hull::Bool=false) where {N}
    vinner = Vector{Vector{N}}(undef, length(dirs))
    underapproximate!(vinner, X, dirs; apply_convex_hull=apply_convex_hull)
end

# in-place version
function underapproximate!(vinner, X::LazySet{N}, dirs::AbstractDirections;
                           apply_convex_hull::Bool=false) where {N}
   @assert dim(X) == dim(dirs) "the dimension of the set, $(dim(X)), doesn't match " *
                               "the dimension of the template directions, $(dim(dirs))"

    @inbounds for (i, di) in enumerate(dirs)
        vinner[i] = σ(di, X)
    end
    if apply_convex_hull
        convex_hull!(vinner)
    end
    return VPolytope(vinner)
end

function underapproximate(X::LazySet, dirs::Type{<:AbstractDirections}; kwargs...)
    return underapproximate(X, dirs(dim(X)), kwargs...)
end
