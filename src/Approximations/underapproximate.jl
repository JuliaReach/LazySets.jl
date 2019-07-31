"""
    underapproximate(X::LazySet{N}, dirs::AbstractDirections) where {N<:Real}

Compute the underapproximation of a convex set by sampling support vectors.

### Input

- `X`    -- set
- `dirs` -- directions

### Output

The `VPolytope` obtained by taking the convex hull of the support vectors of `X`
along the directions determined by `dirs`.
"""
function underapproximate(X::LazySet{N}, dirs::AbstractDirections) where {N<:Real}
    @assert dim(X) == dim(dirs)
    vinner = [σ(di, X) for di in dirs]
    return VPolytope(vinner)
end
