import Base: ∈

export Rectification

"""
    RectificationCache{N<:Real}

Struct that is used as a cache for [`Rectification`](@ref)s.

### Fields

- `set`                -- set represented by the rectification (can be `nothing`
                          if not computed yet)
- `use_support_vector` -- flag indicating whether to use support-vector
                          computations for the cached set
"""
mutable struct RectificationCache{N<:Real}
    set::Union{LazySet{N}, UnionSetArray{N}, Nothing}
    use_support_vector::Bool

    # constructor without a set
    RectificationCache{N}(::Nothing) where {N<:Real} = new{N}(nothing, false)

    # constructor with a set
    RectificationCache{N}(set::LazySet{N}) where {N<:Real} = new{N}(set, true)
end

"""
    Rectification{N<:Real, S<:LazySet{N}}

Type that represents the rectification of a convex set.

### Fields

- `X`     -- convex set
- `cache` -- storage of information computed before

### Notes

Given a vector ``v = (v_1, …, v_n)``, its rectification is defined as
``\\text{rectify}(v) = (v_1', …, v_n')`` such that ``v_i' = \\max(v_i, 0)`` for
each ``i = 1, …, n``.

The extension to a set ``X`` is defined elementwise:

```math
    \\text{rectify}(X) = \\{\\text{rectify}(x) \\mid x ∈ X\\}
```

The rectification of a convex set ``X`` is not necessarily convex.
It can be expressed exactly as the union of the intersection of ``X`` with the
nonnegative orthant and the projection of the intersection of ``X`` with each
other orthant.
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

For example, consider a two-dimensional case and call the orthants
``O_1, …, O_4`` in clockwise fashion, starting with the nonnegative orthant.
We conclude that

```math
    \\text{rectify}(X) = (X ∩ O_1) ∪ \\text{rectify}(X ∩ O_2) ∪ \\text{rectify}(X ∩ O_3) ∪ \\text{rectify}(X ∩ O_4).
```

The rectification of the intersection in the nonpositive orthant,
``\\text{rectify}(X ∩ O_3)``, is either the empty set or the singleton
containing the origin.
The rectification of ``X ∩ O_2`` and ``X ∩ O_4`` both result in flat
``1``-dimensional line segments on the corresponding hyperplane of ``O_1``.
"""
struct Rectification{N<:Real, S<:LazySet{N}}
    X::S
    cache::RectificationCache{N}

    # default constructor that initializes cache
    function Rectification{N, S}(X::S) where {N<:Real, S<:LazySet{N}}
        if X isa AbstractHyperrectangle || X isa CartesianProduct ||
                X isa CartesianProductArray
            # set types with efficient support-vector computations
            set = X
        elseif dim(X) == 1 && isbounded(X)
            # one-dimensional bounded sets are converted to Interval type
            set = convert(Interval, X)
        else
            # do not pre-compute a set; computation is triggered later
            set = nothing
        end
        cache = RectificationCache{N}(set)
        return new{N, S}(X, cache)
    end
end

isoperation(::Type{Rectification}) = true

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

The support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.
"""
function σ(d::AbstractVector{N}, r::Rectification{N}) where {N<:Real}
    if r.cache.set == nothing
        r.cache.set = to_union_of_projections(r)
    end
    return σ(d, r.cache.set)
end

"""
    σ(d::AbstractVector{N},
      r::Rectification{N, <:AbstractHyperrectangle{N}}) where {N<:Real}

Return the support vector of the rectification of a hyperrectangular set.

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

Return the support vector of the rectification of a Cartesian product of two
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

Return the support vector of the rectification of a Cartesian product of a
finite number of convex sets.

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
    ρ(d::AbstractVector{N}, r::Rectification{N}) where {N<:Real}

Evaluate the support function of a rectification of a convex set in a given
direction.

### Input

- `d` -- direction
- `r` -- rectification of a convex set

### Output

The support value of the rectification of a convex set in the given direction.

### Algorithm

We use different procedures for different types of input sets.
If the wrapped set has a suitable structure for which we can efficiently compute
the support vector, we fall back to the evaluation of the support function by
means of the support vector.
Otherwise we compute the union of projections to obtain a precise result (see
[`to_union_of_projections`](@ref)), and then compute the support function for
this union.
(The union is cached internally, so subsequent queries are more efficient.)
"""
function ρ(d::AbstractVector{N}, r::Rectification{N}) where {N<:Real}
    if r.cache.use_support_vector
        return dot(d, σ(d, r))
    end
    if r.cache.set == nothing
        r.cache.set = to_union_of_projections(r)
    end
    return ρ(d, r.cache.set)
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
        d = SingleEntryVector(i, n, one(N))
        if ρ(d, r.X) == N(Inf)
            return false
        end
    end
    return true
end

"""
    to_union_of_projections(r::Rectification{N},
                            concrete_intersection::Bool=false
                           ) where {N<:Real}

Compute an equivalent union of projections from a rectification of a convex set.

### Input

- `r`                     -- rectification of a convex set
- `concrete_intersection` -- (optional, default: `false`) option to compute
                             all intersections concretely or lazily

### Algorithm

Let ``X`` be the set wrapped by the rectification ``r``.
We compute a union of sets that represents the rectification of ``X`` precisely.
The sets are lazy projections, potentially of intersections.

We first identify those dimensions where ``X`` is negative, using one
support-function query per dimension, and collect the dimensions in the index
set ``I_\\text{neg}``.
For each element in ``I_\\text{neg}`` we will later apply a projection to zero.

Next we identify those dimensions from ``I_\\text{neg}`` where ``X`` is also
positive, using another support-function query in each dimension, and collect
the dimensions in the index set ``I_\\text{mix}``.
Let us call the remaining dimensions
(``I_\\text{neg} \\setminus I_\\text{mix}``) ``I_\\text{nonpos}``.
For each dimension in ``j ∈ I_\\text{mix}`` we will apply an intersection with
axis-aligned polyhedra.
In particular, we distinguish two cases using half-spaces ``x_j ≤ 0`` and
``x_j ≥ 0``, and then compute all possible combinations to intersect, using one
half-space per dimension ``j ∈ I_\\text{mix}``.

Next we project the intersections in all dimensions from ``i ∈ I_\\text{mix}``
such that we used the half-space ``x_i ≤ 0`` in their computation, and in all
dimensions ``j ∈ I_\\text{nonpos}`` irrespective of the half-space used.

Finally, we take the union of the resulting sets.

### Output

The result can be one of three cases depending on the wrapped set ``X``, namely
* the set ``X`` if ``X`` is contained in the positive quadrant,
* a `LinearMap` (projection) of ``X`` if for each dimension, ``X`` is only
  either positive or negative, or
* a `UnionSetArray` of `LinearMaps` (projections) otherwise.
"""
function to_union_of_projections(r::Rectification{N},
                                 concrete_intersection::Bool=false
                                ) where {N<:Real}
    n = dim(r.X)

    negative_dimensions, mixed_dimensions, mixed_dimensions_enumeration =
        compute_negative_and_mixed_dimensions(r.X, n)

    if isempty(mixed_dimensions)
        # no mixed dimensions: no intersections needed
        if isempty(negative_dimensions)
            # no negative dimensions: set is unaffected by rectification
            result = r.X
        else
            # project negative dimensions to zero
            v = ones(N, n)
            v[negative_dimensions] .= zero(N)
            result = Diagonal(v) * r.X
        end
    else
        # apply intersections with half-spaces in every mixed dimension
        m = length(mixed_dimensions)
        i = m
        projections = Vector{LinearMap{N}}()
        sizehint!(projections, 2^m)
        while true
            # compute next projection
            polyhedron = construct_constraints(mixed_dimensions,
                                               mixed_dimensions_enumeration, n,
                                               N)
            cap = concrete_intersection ?
                intersection(r.X, polyhedron) : r.X ∩ polyhedron
            projection = construct_projection(cap,
                                              negative_dimensions,
                                              mixed_dimensions,
                                              mixed_dimensions_enumeration, n)
            push!(projections, projection)

            if mixed_dimensions_enumeration[i]
                @inbounds while i > 0 && mixed_dimensions_enumeration[i]
                    i -= 1
                end
                if i == 0
                    break
                end
                mixed_dimensions_enumeration[i] = true
                mixed_dimensions_enumeration[i+1:m] .= false
                i = m
            else
                mixed_dimensions_enumeration[i] = true
            end
        end
        result = UnionSetArray(projections)
    end
    return result
end

# =========================
# internal helper functions
# =========================

# identify dimensions where set is negative and mixed (positive/negative) and
# also fill a bit vector for mixed dimensions
function compute_negative_and_mixed_dimensions(X::LazySet{N}, n) where {N<:Real}
    negative_dimensions = Vector{Int}()
    mixed_dimensions = Vector{Int}()
    mixed_dimensions_enumeration = BitArray(undef, 0)
    for i in 1:n
        d = SingleEntryVector(i, n, -one(N))
        if ρ(d, X) > zero(N)
            push!(negative_dimensions, i)
            if ρ(-d, X) > zero(N)
                push!(mixed_dimensions_enumeration, false)
                push!(mixed_dimensions, i)
            end
        end
    end
    return (negative_dimensions, mixed_dimensions, mixed_dimensions_enumeration)
end

# construct a polyhedron corresponding to an index vector and a bit vector
# (see `to_union_of_projections`)
function construct_constraints(indices, bits, n, ::Type{N}) where {N<:Real}
    P = HPolyhedron{N}()
    for (i, b) in enumerate(bits)
        direction = SingleEntryVector(indices[i], n,
                                                     (b ? -one(N) : one(N)))
        hs = HalfSpace(direction, zero(N))
        addconstraint!(P, hs)
    end
    return P
end

# project negative dimensions to zero (see `to_union_of_projections`)
function construct_projection(set::LazySet{N}, negative_dimensions,
                              mixed_dimensions, mixed_dimensions_enumeration, n
                             ) where {N<:Real}
    # initialize negative indices to zero and all other dimensions to one
    v = ones(N, n)
    v[negative_dimensions] .= zero(N)

    for (i, d) in enumerate(mixed_dimensions)
        if mixed_dimensions_enumeration[i]
            # reset mixed index to one again because the set is nonnegative
            v[d] = one(N)
        end
    end
    return Diagonal(v) * set
end
