export CachedMinkowskiSumArray,
       forget_sets!

# =============================================================
# Minkowski sum of an array of sets with a support vector cache
# =============================================================

"""
    CachedPair{N}

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
    CachedMinkowskiSumArray{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the Minkowski sum of a finite number of convex sets.
Support vector queries are cached.

### Fields

- `array` -- array of convex sets
- `cache` -- cache of support vector query results

### Notes

This type assumes that the dimensions of all elements match.

The `ZeroSet` is the neutral element and the `EmptySet` is the absorbing element
for `CachedMinkowskiSumArray`.

The cache (field `cache`) is implemented as dictionary whose keys are directions
and whose values are pairs `(k, s)` where `k` is the number of elements in the
array `array` when the support vector was evaluated last time, and `s` is the
support vector that was obtained.
Thus this type assumes that `array` is not modified except by adding new sets at
the end.

Constructors:

- `CachedMinkowskiSumArray(array::Vector{<:LazySet})` -- default constructor

- `CachedMinkowskiSumArray([n]::Int=0, [N]::Type=Float64)`
  -- constructor for an empty sum with optional size hint and numeric type
"""
struct CachedMinkowskiSumArray{N<:Real, S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
    cache::Dict{AbstractVector{N}, CachedPair{N}}

    # default constructor that initializes cache
    CachedMinkowskiSumArray{N, S}(arr::Vector{S}) where {N<:Real, S<:LazySet{N}} =
        new{N, S}(arr, Dict{AbstractVector{N}, CachedPair{N}}())
end

isoperationtype(::Type{<:CachedMinkowskiSumArray}) = true
isconvextype(::Type{CachedMinkowskiSumArray{N, S}}) where {N, S} = isconvextype(S)

# convenience constructor without type parameter
CachedMinkowskiSumArray(arr::Vector{S}) where {N<:Real, S<:LazySet{N}} =
    CachedMinkowskiSumArray{N, S}(arr)

# constructor for an empty sum with optional size hint and numeric type
function CachedMinkowskiSumArray(n::Int=0, N::Type=Float64)
    arr = Vector{LazySet{N}}()
    sizehint!(arr, n)
    return CachedMinkowskiSumArray(arr)
end

# ZeroSet is the neutral element for CachedMinkowskiSumArray
@neutral(CachedMinkowskiSumArray, ZeroSet)

# EmptySet and Universe are the absorbing element for CachedMinkowskiSumArray
@absorbing(CachedMinkowskiSumArray, EmptySet)
# @absorbing(CachedMinkowskiSumArray, Universe)  # TODO problematic

"""
    array(cms::CachedMinkowskiSumArray{N, S}) where {N<:Real, S<:LazySet{N}}

Return the array of a caching Minkowski sum.

### Input

- `cms` -- caching Minkowski sum

### Output

The array of a caching Minkowski sum.
"""
function array(cms::CachedMinkowskiSumArray{N, S}
              ) where {N<:Real, S<:LazySet{N}}
    return cms.array
end

"""
    dim(cms::CachedMinkowskiSumArray)

Return the dimension of a caching Minkowski sum.

### Input

- `cms` -- caching Minkowski sum

### Output

The ambient dimension of the caching Minkowski sum.
"""
function dim(cms::CachedMinkowskiSumArray)
    return length(cms.array) == 0 ? 0 : dim(cms.array[1])
end

"""
    σ(d::AbstractVector{N}, cms::CachedMinkowskiSumArray{N}) where {N<:Real}

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
function σ(d::AbstractVector{N}, cms::CachedMinkowskiSumArray{N}) where {N<:Real}
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
            svec = svec1 + σ_helper(d, @view arr[k+1:l])
        end
    else
        # first-time computation of support vector
        svec = σ_helper(d, arr)
    end
    # NOTE: make a copy of the direction vector (can be modified outside)
    cache[copy(d)] = CachedPair(l, svec)
    return svec
end

"""
	isbounded(cms::CachedMinkowskiSumArray)

Determine whether a caching Minkowski sum is bounded.

### Input

- `cms` -- caching Minkowski sum

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(cms::CachedMinkowskiSumArray)
    return all(x -> isbounded(x), cms.array)
end

"""
    isempty(cms::CachedMinkowskiSumArray)

Return if a caching Minkowski sum array is empty or not.

### Input

- `cms` -- caching Minkowski sum

### Output

`true` iff any of the wrapped sets are empty.

### Notes

Forgotten sets cannot be checked anymore.
Usually they have been empty because otherwise the support vector query should
have crashed before.
In that case, the caching Minkowski sum should not be used further.
"""
function isempty(cms::CachedMinkowskiSumArray)
    return any(X -> isempty(X), array(cms))
end

"""
    forget_sets!(cms::CachedMinkowskiSumArray)

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

julia> cms1 = CachedMinkowskiSumArray(2); cms2 = CachedMinkowskiSumArray(2);

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
function forget_sets!(cms::CachedMinkowskiSumArray)
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
