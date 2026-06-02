"""
    σ(d::AbstractVector, cms::CachedMinkowskiSumArray)

Return a support vector of a cached Minkowski sum in a given direction.

### Input

- `d`   -- direction
- `cms` -- cached Minkowski sum

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.

### Notes

The result is cached, i.e., any further query with the same direction runs in
constant time.
When sets are added to the cached Minkowski sum, the query is only performed
for the new sets.
"""
@validate function σ(d::AbstractVector, cms::CachedMinkowskiSumArray)
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
            svec = svec1 + _σ_msum_array(d, @view arr[(k + 1):l])
        end
    else
        # first-time computation of support vector
        svec = _σ_msum_array(d, arr)
    end
    # NOTE: make a copy of the direction vector (can be modified outside)
    cache[copy(d)] = CachedPair(l, svec)
    return svec
end
