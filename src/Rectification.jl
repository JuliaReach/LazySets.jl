export Rectification

"""
    RectificationCache{N<:Real}

Struct that is used as a cache for [`Rectification`](@ref)s.
"""
struct RectificationCache{N<:Real}
    # empty for now
end

"""
    Rectification{N<:Real, S<:LazySet{N}}

Type that represents the rectification of a convex set.

### Fields

- `X`     -- convex set
- `cache` -- storage of information computed before

### Notes

The rectification of a convex set ``X`` is not necessarily convex.
It can be expressed exactly as the union of the intersection of ``X`` with the
nonnegative orthant and the projection of the intersection of ``X`` with each
orthant where exactly one dimension is nonpositive on the corresponding
hyperplane of the nonnegative orthant.
This can be seen as follows.

First we observe that rectification distributes with union.

```math
    \\text{rectify}(X_1 ∪ … ∪ X_m) = ⋃_j \\text{rectify}(X_j)
```

Next we express ``X`` as the union of the intersection of ``X`` with each
orthant ``O``.

```math
    X = ⋃_j (X ∩ O_j)
```

Thus we have

```math
    \\text{rectify}(X) = \\text{rectify}((X ∩ O_1) ∪ … ∪ (X ∩ O_m)) = ⋃_j \\text{rectify}(X ∩ O_j).
```

Clearly, ``\\text{rectify}(X ∩ O_j) = X`` if ``O_j`` is the nonnegative orthant.

Furthermore, we need not consider orthants with two or more negative dimensions
because their contributions to the final set are subsumed.

For example, consider a two-dimensional case and call the orthants
``O_1, …, O_4`` in clockwise fashion, starting with the nonnegative orthant.
The rectification of the intersection in the nonpositive orthant,
``\\text{rectify}(X ∩ O_3)``, is either the empty set or the singleton
containing the origin; in the latter case, at least one of
``\\text{rectify}(X ∩ O_2)`` and ``\\text{rectify}(X ∩ O_4)`` also contain the
origin.
We conclude that

```math
    \\text{rectify}(X) = (X ∩ O_1) ∪ \\text{rectify}(X ∩ O_2) ∩ \\text{rectify}(X ∩ O_4)
```

and that all applications of ``\\text{rectify}`` on the right-hand side result
in flat ``n-1``-dimensional shapes on the corresponding hyperplane of ``O_1``.
"""
struct Rectification{N<:Real, S<:LazySet{N}}
    X::S
    cache::RectificationCache{N}

    # default constructor that initializes cache
    function Rectification{N, S}(X::S) where {N<:Real, S<:LazySet{N}}
        cache = RectificationCache{N}()
        return new{N, S}(X, cache)
    end
end

# convenience constructor without type parameter
Rectification(X::S) where {N<:Real, S<:LazySet{N}} = Rectification{N, S}(X)
