"""
    RectificationCache{N}

Struct that is used as a cache for [`Rectification`](@ref)s.

### Fields

- `set`                -- set represented by the rectification (can be `nothing`
                          if not computed yet)
- `use_support_vector` -- flag indicating whether to use support-vector
                          computations for the cached set
"""
mutable struct RectificationCache{N}
    set::Union{LazySet{N},Nothing}
    use_support_vector::Bool

    # constructor without a set
    RectificationCache{N}(::Nothing) where {N} = new{N}(nothing, false)

    # constructor with a set
    RectificationCache{N}(set::LazySet{N}) where {N} = new{N}(set, true)
end

"""
    Rectification{N, S<:LazySet{N}} <: LazySet{N}

Type that represents the rectification of a set.

### Fields

- `X`     -- set
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
struct Rectification{N,S<:LazySet{N}} <: LazySet{N}
    X::S
    cache::RectificationCache{N}

    # default constructor that initializes cache
    function Rectification(X::S) where {N,S<:LazySet{N}}
        if X isa AbstractHyperrectangle || X isa CartesianProduct ||
           X isa CartesianProductArray || X isa EmptySet
            # set types with efficient support-vector computations
            set = X
        elseif isconvextype(S) && dim(X) == 1 && isbounded(X)
            # one-dimensional compact convex sets are converted to Interval type
            set = convert(Interval, X)
        else
            # do not pre-compute a set; computation is triggered later
            set = nothing
        end
        cache = RectificationCache{N}(set)
        return new{N,S}(X, cache)
    end
end

function _compute_exact_representation!(R::Rectification)
    if isnothing(R.cache.set)
        R.cache.set = to_union_of_projections(R)
    end
    return R.cache.set
end

"""
    set(R::Rectification)

Return the original set of a rectification.

### Input

- `R` -- rectification

### Output

The original set of the rectification.
"""
function set(R::Rectification)
    return R.X
end

"""
    to_union_of_projections(R::Rectification,
                            [concrete_intersection]::Bool=false) where

Compute an equivalent union of projections from a rectification.

### Input

- `R`                     -- rectification
- `concrete_intersection` -- (optional, default: `false`) option to compute
                             all intersections concretely or lazily
- `filter_empty_sets`     -- (optional, default: `true`) option to filter out
                             empty sets in the union

### Algorithm

Let ``X`` be the set wrapped by the rectification ``R``.
We compute a union of sets that represents the rectification of ``X`` precisely.
The sets are lazy projections, potentially of intersections.

We first identify those dimensions where ``X`` is negative, using one call to
`low` per dimension, and collect the dimensions in the index set
``I_\\text{neg}``.
For each element in ``I_\\text{neg}`` we will later apply a projection to zero.

Next we identify those dimensions from ``I_\\text{neg}`` where ``X`` is also
positive, using another `high` query in each dimension, and collect the
dimensions in the index set ``I_\\text{mix}``.
Let us call the remaining dimensions
(``I_\\text{neg} ∖ I_\\text{mix}``) ``I_\\text{nonpos}``.
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
* a `UnionSetArray` of `LinearMap`s (projections) otherwise.
"""
function to_union_of_projections(R::Rectification,
                                 concrete_intersection::Bool=false;
                                 filter_empty_sets::Bool=true)
    N = eltype(R)
    n = dim(R.X)

    negative_dimensions, mixed_dimensions, mixed_dimensions_enumeration = compute_negative_and_mixed_dimensions(R.X,
                                                                                                                n)

    if isempty(mixed_dimensions)
        # no mixed dimensions: no intersections needed
        if isempty(negative_dimensions)
            # no negative dimensions: set is unaffected by rectification
            result = R.X
        else
            # project negative dimensions to zero
            v = ones(N, n)
            v[negative_dimensions] .= zero(N)
            result = Diagonal(v) * R.X
        end
    else
        # apply intersections with half-spaces in every mixed dimension
        m = length(mixed_dimensions)
        i = m
        projections = Vector{LinearMap{N}}()
        while true
            # compute next projection
            polyhedron = construct_constraints(mixed_dimensions,
                                               mixed_dimensions_enumeration, n,
                                               N)
            if !filter_empty_sets || !isdisjoint(R.X, polyhedron)
                if concrete_intersection
                    cap = intersection(R.X, polyhedron)
                else
                    cap = Intersection(R.X, polyhedron)
                    if filter_empty_sets
                        set_isempty!(cap, false)  # intersection is known to not be empty
                    end
                end
                projection = construct_projection(cap,
                                                  negative_dimensions,
                                                  mixed_dimensions,
                                                  mixed_dimensions_enumeration, n)
                push!(projections, projection)
            end

            if mixed_dimensions_enumeration[i]
                @inbounds while i > 0 && mixed_dimensions_enumeration[i]
                    i -= 1
                end
                if i == 0
                    break
                end
                mixed_dimensions_enumeration[i] = true
                mixed_dimensions_enumeration[(i + 1):m] .= false
                i = m
            else
                mixed_dimensions_enumeration[i] = true
            end
        end
        result = UnionSetArray(projections)
    end
    return result
end

# identify dimensions where the set is negative and mixed (positive/negative)
# and also fill a bit vector for mixed dimensions
function compute_negative_and_mixed_dimensions(X::LazySet{N}, n) where {N}
    negative_dimensions = Vector{Int}()
    mixed_dimensions = Vector{Int}()
    mixed_dimensions_enumeration = BitArray(undef, 0)
    for i in 1:n
        if low(X, i) < zero(N)
            push!(negative_dimensions, i)
            if high(X, i) > zero(N)
                push!(mixed_dimensions_enumeration, false)
                push!(mixed_dimensions, i)
            end
        end
    end
    return (negative_dimensions, mixed_dimensions, mixed_dimensions_enumeration)
end

# construct a polyhedron corresponding to an index vector and a bit vector
# (see `to_union_of_projections`)
function construct_constraints(indices, bits, n, ::Type{N}) where {N}
    P = HPolyhedron{N,SingleEntryVector{N}}()
    for (i, b) in enumerate(bits)
        a = SingleEntryVector(indices[i], n, (b ? -one(N) : one(N)))
        hs = HalfSpace(a, zero(N))
        addconstraint!(P, hs)
    end
    return P
end

# project negative dimensions to zero (see `to_union_of_projections`)
function construct_projection(X::LazySet{N}, negative_dimensions,
                              mixed_dimensions, mixed_dimensions_enumeration, n) where {N}
    # initialize negative indices to zero and all other dimensions to one
    v = ones(N, n)
    v[negative_dimensions] .= zero(N)

    for (i, d) in enumerate(mixed_dimensions)
        if mixed_dimensions_enumeration[i]
            # reset mixed index to one again because the set is nonnegative
            v[d] = one(N)
        end
    end
    return Diagonal(v) * X
end

include("an_element.jl")
include("concretize.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isconvextype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("in.jl")
include("support_function.jl")
include("support_vector.jl")
