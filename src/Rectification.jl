import Base: ∈

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

"""
    dim(r::Rectification)::Int

Return the dimension of a rectification.

### Input

- `r` -- rectification

### Output

The ambient dimension of the rectification.
"""
function dim(r::Rectification)::Int
    return dim(r.X)
end

"""
    σ(d::AbstractVector{N}, r::Rectification{N}) where {N<:Real}

Return the support vector of a rectification.

### Input

- `d` -- direction
- `r` -- rectification

### Output

If the rectified set is one-dimensional, we convert it to an `Interval` and call
a specialized method.
Otherise, we throw an error because we currently do not support higher
dimensions.
"""
function σ(d::AbstractVector{N}, r::Rectification{N}) where {N<:Real}
    if dim(r) == 1
        # rectification in 1D is easy
        r2 = Rectification(convert(Interval, r.X))
        return σ(d, r2)
    end
    error("the exact support vector of a rectification is not implemented")
end

"""
    σ(d::AbstractVector{N},
      r::Rectification{N, <:AbstractHyperrectangle{N}}) where {N<:Real}

Return the support vector of a rectification of a hyperrectangular set.

### Input

- `d` -- direction
- `r` -- rectification of a hyperrectangular set

### Output

The support vector in the given direction.

### Algorithm

Let ``r(·)`` be the rectification of a vector respectively a set, and let ``H``
be a hyperrectangle.
Then ``σ_{r(H)}(d) = r(σ_{H}(d))``.
"""
function σ(d::AbstractVector{N},
           r::Rectification{N, <:AbstractHyperrectangle{N}}) where {N<:Real}
    return rectify(σ(d, r.X))
end

"""
    σ(d::AbstractVector{N},
      r::Rectification{N, <:CartesianProduct{N}}) where {N<:Real}

Return the support vector of a rectification of a Cartesian product of two
convex sets.

### Input

- `d` -- direction
- `r` -- rectification of a Cartesian product of two convex sets

### Output

The support vector in the given direction.

### Algorithm

Rectification distributes with the Cartesian product.
Let ``r(·)`` be the rectification of a set.
We can just query the support vector for ``r(X)`` and ``r(Y)`` recursively:
``σ_{r(X × Y)}(d) = σ_{r(X)}(d_X) × σ_{r(Y)}(d_Y)``, where ``x × y``
concatenates vectors ``x`` and ``y``.
"""
function σ(d::AbstractVector{N},
           r::Rectification{N, <:CartesianProduct{N}}) where {N<:Real}
    X, Y = r.X.X, r.X.Y
    n1 = dim(X)
    return vcat(σ(d[1:n1], Rectification(X)), σ(d[n1+1:end], Rectification(Y)))
end

"""
    σ(d::AbstractVector{N},
      r::Rectification{N, <:CartesianProductArray{N}}) where {N<:Real}

Return the support vector of a rectification of a Cartesian product of a finite
number of convex sets.

### Input

- `d` -- direction
- `r` -- rectification of a Cartesian product of a finite number of convex sets

### Output

The support vector in the given direction.

### Algorithm

Rectification distributes with the Cartesian product.
Let ``r(·)`` be the rectification of a set.
We can just query the support vector for each subspace recursively:
``σ_{r(X_1 × ⋯ × X_m)}(d) = σ_{r(X_1)}(d_{X_1}) × ⋯ × σ_{r(X_m)}(d_{X_m})``,
where ``x × y`` concatenates vectors ``x`` and ``y``.
"""
function σ(d::AbstractVector{N},
           r::Rectification{N, <:CartesianProductArray{N}}) where {N<:Real}
    svec = similar(d)
    i = 1
    for X in array(r.X)
        nX = dim(X)
        j = i + nX - 1
        svec[i:j] = σ(d[i:j], Rectification(X))
        i = j + 1
    end
    return svec
end

"""
    an_element(r::Rectification{N})::Vector{N} where {N<:Real}

Return some element of a rectification.

### Input

- `r` -- rectification

### Output

An element in the rectification.
The implementation relies on the `an_element` function of the wrapped set.
"""
function an_element(r::Rectification{N})::Vector{N} where {N<:Real}
    return rectify(an_element(r.X))
end

"""
    ∈(x::AbstractVector{N}, r::Rectification{N})::Bool where {N<:Real}

Check whether a given point is contained in a rectification.

### Input

- `x` -- point/vector
- `r` -- rectification

### Output

`true` iff ``x ∈ r``.

### Algorithm

We first scan for negative entries in the vector.
If there are any, the vector is not contained in the rectification.

Next we ask a membership query in the wrapped set.
If the answer is positive, the vector is contained in the rectification.

Otherwise, we scan for zero entries in the vector.
If there are none, membership reduces to membership in the wrapped set, and so
the answer is negative.

Finally, if there are zero entries in the vector and the vector is not contained
in the wrapped set, we give up and throw an error.
"""
function ∈(x::AbstractVector{N}, r::Rectification{N})::Bool where {N<:Real}
    # scan for negative entries
    if any(xi -> xi < zero(N), x)
        return false
    end

    # membership test in the wrapped set
    if x ∈ r.X
        return true
    end

    # scan for zero entries
    if all(xi -> !iszero(xi), x)
        return false
    end

    # give up
    error("cannot determine membership of a vector with zero entries in a " *
          "lazy rectification")
end

"""
    isempty(r::Rectification)::Bool

Check whether a rectification is empty or not.

### Input

- `r` -- rectification

### Output

`true` iff the wrapped set is empty.
"""
function isempty(r::Rectification)::Bool
    return isempty(r.X)
end

"""
    isbounded(r::Rectification)::Bool

Determine whether a rectification is bounded.

### Input

- `r` -- rectification

### Output

`true` iff the rectification is bounded.

### Algorithm

Let ``X`` be the set wrapped by rectification ``r``.
We first check whether ``X`` is bounded (because then ``r`` is bounded).
Otherwise, we check unboundedness of ``X`` in direction ``(1, 1, …, 1)``, which
is sufficient for unboundedness of ``r``; this step is not necessary but rather
a heuristics.
Otherwise, we check boundedness of ``X`` in every positive unit direction, which
is sufficient and necessary for boundedness of ``r``.
"""
function isbounded(r::Rectification{N})::Bool where {N<:Real}
    # check boundedness of the wrapped set
    if isbounded(r.X)
        return true
    end

    # heuristics: check boundedness in the one-vector direction
    n = dim(r.X)
    if ρ(ones(N, n), r.X) == N(Inf)
        return false
    end

    # check boundedness in positive unit directions
    @inbounds for i in 1:n
        d = Approximations.UnitVector(i, n, one(N))
        if ρ(d, r.X) == N(Inf)
            return false
        end
    end
    return true
end
