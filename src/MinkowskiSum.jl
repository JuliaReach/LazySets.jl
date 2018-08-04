import Base:+, getindex

export MinkowskiSum, ⊕,
       MinkowskiSumArray,
       MinkowskiSum!,
       array,
       CacheMinkowskiSum,
       forget_sets!

"""
    MinkowskiSum{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the Minkowski sum of two convex sets.

### Fields

- `X` -- first convex set
- `Y` -- second convex set

### Notes

The `ZeroSet` is the neutral element and the `EmptySet` is the absorbing element
for `MinkowskiSum`.
"""
struct MinkowskiSum{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function MinkowskiSum{N, S1, S2}(X::S1, Y::S2) where
            {N<:Real, S1<:LazySet{N}, S2<:LazySet{N}}
        @assert dim(X) == dim(Y) "sets in a Minkowski sum must have the " *
            "same dimension"
        return new(X, Y)
    end
end
# type-less convenience constructor
MinkowskiSum(X::S1, Y::S2) where {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
    MinkowskiSum{N, S1, S2}(X, Y)

# ZeroSet is the neutral element for MinkowskiSum
@neutral(MinkowskiSum, ZeroSet)

# EmptySet is the absorbing element for MinkowskiSum
@absorbing(MinkowskiSum, EmptySet)

"""
    X + Y

Convenience constructor for Minkowski sum.

### Input

- `X` -- a convex set
- `Y` -- another convex set

### Output

The symbolic Minkowski sum of ``X`` and ``Y``.
"""
+(X::LazySet, Y::LazySet) = MinkowskiSum(X, Y)

"""
    ⊕(X::LazySet, Y::LazySet)

Unicode alias constructor ⊕ (`oplus`) for the lazy Minkowski sum operator.
"""
⊕(X::LazySet, Y::LazySet) = MinkowskiSum(X, Y)

"""
    dim(ms::MinkowskiSum)::Int

Return the dimension of a Minkowski sum.

### Input

- `ms` -- Minkowski sum

### Output

The ambient dimension of the Minkowski sum.
"""
function dim(ms::MinkowskiSum)::Int
    return dim(ms.X)
end

"""
    σ(d::AbstractVector{N}, ms::MinkowskiSum{N}) where {N<:Real}

Return the support vector of a Minkowski sum.

### Input

- `d`  -- direction
- `ms` -- Minkowski sum

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.
"""
function σ(d::AbstractVector{N}, ms::MinkowskiSum{N}) where {N<:Real}
    return σ(d, ms.X) + σ(d, ms.Y)
end

# =================================
# Minkowski sum of an array of sets
# =================================

"""
    MinkowskiSumArray{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the Minkowski sum of a finite number of convex sets.

### Fields

- `array` -- array of convex sets

### Notes

This type assumes that the dimensions of all elements match.

The `ZeroSet` is the neutral element and the `EmptySet` is the absorbing element
for `MinkowskiSumArray`.

Constructors:

- `MinkowskiSumArray(array::Vector{<:LazySet})` -- default constructor

- `MinkowskiSumArray([n]::Int=0, [N]::Type=Float64)`
  -- constructor for an empty sum with optional size hint and numeric type
"""
struct MinkowskiSumArray{N<:Real, S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

# type-less convenience constructor
MinkowskiSumArray(arr::Vector{S}) where {S<:LazySet{N}} where {N<:Real} =
    MinkowskiSumArray{N, S}(arr)

# constructor for an empty sum with optional size hint and numeric type
function MinkowskiSumArray(n::Int=0, N::Type=Float64)::MinkowskiSumArray
    arr = Vector{LazySet{N}}(0)
    sizehint!(arr, n)
    return MinkowskiSumArray(arr)
end

# ZeroSet is the neutral element for MinkowskiSumArray
@neutral(MinkowskiSumArray, ZeroSet)

# EmptySet is the absorbing element for MinkowskiSumArray
@absorbing(MinkowskiSumArray, EmptySet)

# add functions connecting MinkowskiSum and MinkowskiSumArray
@declare_array_version(MinkowskiSum, MinkowskiSumArray)

"""
    array(msa::MinkowskiSumArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}

Return the array of a Minkowski sum of a finite number of convex sets.

### Input

- `msa` -- Minkowski sum array

### Output

The array of a Minkowski sum of a finite number of convex sets.
"""
function array(msa::MinkowskiSumArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}
    return msa.array
end

"""
    dim(msa::MinkowskiSumArray)::Int

Return the dimension of a Minkowski sum of a finite number of sets.

### Input

- `msa` -- Minkowski sum array

### Output

The ambient dimension of the Minkowski sum of a finite number of sets.
"""
function dim(msa::MinkowskiSumArray)::Int
    return length(msa.array) == 0 ? 0 : dim(msa.array[1])
end

"""
    σ(d::AbstractVector{N}, msa::MinkowskiSumArray{N}) where {N<:Real}

Return the support vector of a Minkowski sum of a finite number of sets in a
given direction.

### Input

- `d`   -- direction
- `msa` -- Minkowski sum array

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.
"""
function σ(d::AbstractVector{N}, msa::MinkowskiSumArray{N}) where {N<:Real}
    return _σ_helper(d, msa.array)
end

# =============================================================
# Minkowski sum of an array of sets with a support vector cache
# =============================================================

"""
    CachedPair{N} where N

A mutable pair of an index and a vector.

### Fields

- `idx` -- index
- `vec` -- vector
"""
mutable struct CachedPair{N}
    idx::Int
    vec::AbstractVector{N}
end

function getindex(cp::CachedPair, idx::Int)
    if idx == 1
        return cp.idx
    elseif idx == 2
        return cp.vec
    end
    error("invalid index $idx, can only access a pair at index 1 or 2")
end

"""
    CacheMinkowskiSum{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the Minkowski sum of a finite number of convex sets.
Support vector queries are cached.

### Fields

- `array` -- array of convex sets
- `cache` -- cache of support vector query results

### Notes

This type assumes that the dimensions of all elements match.

The `ZeroSet` is the neutral element and the `EmptySet` is the absorbing element
for `CacheMinkowskiSum`.

The cache (field `cache`) is implemented as dictionary whose keys are directions
and whose values are pairs `(k, s)` where `k` is the number of elements in the
array `array` when the support vector was evaluated last time, and `s` is the
support vector that was obtained.
Thus this type assumes that `array` is not modified except by adding new sets at
the end.

Constructors:

- `CacheMinkowskiSum(array::Vector{<:LazySet})` -- default constructor

- `CacheMinkowskiSum([n]::Int=0, [N]::Type=Float64)`
  -- constructor for an empty sum with optional size hint and numeric type
"""
struct CacheMinkowskiSum{N<:Real, S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
    cache::Dict{AbstractVector{N}, CachedPair{N}}

    # default constructor that initializes cache
    CacheMinkowskiSum{N, S}(arr::Vector{S}) where {N<:Real, S<:LazySet{N}} =
        new{N, S}(arr, Dict{AbstractVector{N}, CachedPair{N}}())
end

# type-less convenience constructor
CacheMinkowskiSum(arr::Vector{S}) where {N<:Real, S<:LazySet{N}} =
    CacheMinkowskiSum{N, S}(arr)

# constructor for an empty sum with optional size hint and numeric type
function CacheMinkowskiSum(n::Int=0, N::Type=Float64)::CacheMinkowskiSum
    arr = Vector{LazySet{N}}(0)
    sizehint!(arr, n)
    return CacheMinkowskiSum(arr)
end

# ZeroSet is the neutral element for CacheMinkowskiSum
@neutral(CacheMinkowskiSum, ZeroSet)

# EmptySet is the absorbing element for CacheMinkowskiSum
@absorbing(CacheMinkowskiSum, EmptySet)

"""
    array(cms::CacheMinkowskiSum{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}

Return the array of a caching Minkowski sum.

### Input

- `cms` -- caching Minkowski sum

### Output

The array of a caching Minkowski sum.
"""
function array(cms::CacheMinkowskiSum{N, S}
              )::Vector{S} where {N<:Real, S<:LazySet{N}}
    return cms.array
end

"""
    dim(cms::CacheMinkowskiSum)::Int

Return the dimension of a caching Minkowski sum.

### Input

- `cms` -- caching Minkowski sum

### Output

The ambient dimension of the caching Minkowski sum.
"""
function dim(cms::CacheMinkowskiSum)::Int
    return length(cms.array) == 0 ? 0 : dim(cms.array[1])
end

"""
    σ(d::AbstractVector{N}, cms::CacheMinkowskiSum{N}) where {N<:Real}

Return the support vector of a caching Minkowski sum in a given direction.

### Input

- `d`   -- direction
- `cms` -- caching Minkowski sum

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.

### Notes

The result is cached, i.e., any further query with the same direction runs in
constant time.
When sets are added to the caching Minkowski sum, the query is only performed
for the new sets.
"""
function σ(d::AbstractVector{N}, cms::CacheMinkowskiSum{N}) where {N<:Real}
    arr = array(cms)
    l = length(arr)
    cache = cms.cache
    if haskey(cache, d)
        pair = cache[d]
        k = pair.idx
        svec1 = pair.vec
        if k == l
            # has already stored the full support vector
            return svec1
        else
            # has only stored the support vector of the first k sets
            @assert k < l "invalid cache index"
            svec = svec1 + _σ_helper(d, @view arr[k+1:l])
        end
    else
        # first-time computation of support vector
        svec = _σ_helper(d, arr)
    end
    # NOTE: make a copy of the direction vector (can be modified outside)
    cache[copy(d)] = CachedPair(l, svec)
    return svec
end

"""
    forget_sets!(cms::CacheMinkowskiSum)::Int

Tell a caching Minkowski sum to forget the stored sets (but not the support
vectors).
Only those sets are forgotten such that for each cached direction the support
vector has been computed before.

### Input

- `cms` -- caching Minkowski sum

### Output

The number of sets that have been forgotten.

### Notes

This function should only be used under the assertion that no new directions are
queried in the future; otherwise such support vector results will be incorrect.

This implementation is optimistic and first tries to remove all sets.
However, it also checks that for all cached directions the support vector has
been computed before.
If it finds that this is not the case, the implementation identifies the
biggest index ``k`` such that the above holds for the ``k`` oldest sets, and
then it only removes these.
See the example below.

### Examples

```jldoctest
julia> x1 = BallInf(ones(3), 3.); x2 = Ball1(ones(3), 5.);

julia> cms1 = CacheMinkowskiSum(2); cms2 = CacheMinkowskiSum(2);

julia> d = ones(3);

julia> a1 = array(cms1); a2 = array(cms2);

julia> push!(a1, x1); push!(a2, x1);

julia> σ(d, cms1); σ(d, cms2);

julia> push!(a1, x2); push!(a2, x2);

julia> σ(d, cms1);

julia> idx1 = forget_sets!(cms1) # support vector was computed for both sets
2

julia> idx1 = forget_sets!(cms2) # support vector was only computed for first set
1
```
"""
function forget_sets!(cms::CacheMinkowskiSum)::Int
    len = length(cms.array)
    sets_to_remove = len
    for pair in values(cms.cache)
        if pair.idx != len
            @assert pair.idx < len "invalid cache index"
            sets_to_remove = min(sets_to_remove, pair.idx)
            pair.idx -= len
        else
            pair.idx = 0
        end
    end

    if sets_to_remove != len
        # not all sets shall be forgotten
        offset = len - sets_to_remove
        for pair in values(cms.cache)
            # fix the wrong index assignments from above
            pair.idx += offset
        end
        # forget only the first k sets
        splice!(cms.array, 1:sets_to_remove)
    else
        # forget all sets
        empty!(cms.array)
    end
    return sets_to_remove
end

# ================
# Helper functions
# ================

@inline function _σ_helper(d::AbstractVector{<:Real},
                           array::AbstractVector{<:LazySet})::Vector{<:Real}
    svec = zeros(eltype(d), length(d))
    for sj in array
        svec += σ(d, sj)
    end
    return svec
end
