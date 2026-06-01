"""
    dim(cms::CachedMinkowskiSumArray)

Return the dimension of a cached Minkowski sum.

### Input

- `cms` -- cached Minkowski sum

### Output

The ambient dimension of the cached Minkowski sum, or `0` if there is no set in
the array.
"""
function dim(cms::CachedMinkowskiSumArray)
    return length(cms.array) == 0 ? 0 : dim(cms.array[1])
end
