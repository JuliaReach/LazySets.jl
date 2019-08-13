export hausdorff_distance

"""
    hausdorff_distance(X::LazySet{N}, Y::LazySet{N}; [p]::N=N(Inf),
                       [ε]::N=N(1e-3)) where {N<:Real}

Compute the Hausdorff distance between two convex sets up to a given threshold.

### Input

- `X` -- convex set
- `Y` -- convex set
- `p` -- (optional, default: `Inf`) norm parameter of the Hausdorff distance
- `ε` -- (optional, default: `1e-3`) precision threshold; the true Hausdorff
         distance may diverge from the result by at most this value

### Output

A value from the ``ε``-neighborhood of the Hausdorff distance between ``X`` and
``Y``.

### Notes

Given a ``p``-norm, the Hausdorff distance ``d_H^p(X, Y)`` between sets ``X``
and ``Y`` is defined as follows:

```math
    d_H^p(X, Y) = \\inf\\{δ ≥ 0 \\mid Y ⊆ X ⊕ δ 𝐵_p^n \\text{ and } X ⊆ Y ⊕ δ 𝐵_p^n\\}
```

Here ``𝐵_p^n`` is the ``n``-dimensional unit ball in the ``p``-norm.

The implementation may internally rely on the support function of ``X`` and
``Y``; hence any imprecision in the implementation of the support function may
affect the result.
At the time of writing, the only set type with imprecise support function is the
lazy [`Intersection`](@ref).

### Algorithm

We perform binary search for bounding the Hausdorff distance in an interval
``[l, u]``, where initially ``l`` is ``0`` and ``u`` is described below.
The binary search terminates when ``u - l ≤ ε``, i.e., the interval becomes
sufficiently small.

To find an upper bound ``u``, we start with the heuristics of taking the biggest
distance in the axis-parallel directions.
As long as this bound does not work, we increase the bound by ``2``.

Given a value ``δ``, to check whether the sets are within Hausdorff distance
``δ``, we simply check the inclusions given above, where on the right-hand side
we use a lazy `MinkowskiSum` with a `Ballp` centered in the origin.
"""
function hausdorff_distance(X::LazySet{N}, Y::LazySet{N}; p::N=N(Inf),
                            ε::N=N(1e-3)) where {N<:Real}
    @assert isbounded(X) && isbounded(Y) "the Hausdorff distance is only " *
        "defined for compact sets"

    n = dim(X)
    @assert dim(Y) == n "the Hausdorff distance is only defined between sets " *
                        "of the same dimension, but they had dimensions $n " *
                        "resp. $(dim(Y))"

    # phase 1: find a finite upper bound
    δ_upper = maximum(d -> abs(ρ(d, X) - ρ(d, Y)), BoxDirections{N}(n))
    # verify that this is an upper bound
    while !_mutual_issubset_in_δ_bloating(X, Y, δ_upper, n, p)
        δ_upper *= N(2)
    end

    # phase 2: perform binary search between lower bound (initially 0) and upper
    # bound until convergence
    δ_lower = N(0)
    while δ_upper - δ_lower > ε
        δ = (δ_upper + δ_lower) / N(2)
        if _mutual_issubset_in_δ_bloating(X, Y, δ, n, p)
            δ_upper = δ
        else
            δ_lower = δ
        end
    end
    return δ_upper
end

function _mutual_issubset_in_δ_bloating(X, Y, δ, n, p)
    return _issubset_in_δ_bloating(X, Y, δ, n, p) &&
           _issubset_in_δ_bloating(Y, X, δ, n, p)
end

function _issubset_in_δ_bloating(X::LazySet{N}, Y, δ, n, p) where {N}
    return X ⊆ Y + Ballp(p, zeros(N, n), δ)
end

# for polytopes the default implementation of `⊆` requires membership in the rhs
# set, which will be a MinkowskiSum and hence not available; we use the
# alternative based on constraints_list on the right instead
function _issubset_in_δ_bloating(X::AbstractPolytope{N}, Y, δ, n, p
                                ) where {N<:Real}
    return LazySets._issubset_constraints_list(X, Y + Ballp(p, zeros(N, n), δ))
end


