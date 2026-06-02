"""
    isempty(cms::CachedMinkowskiSumArray)

Check whether a cached Minkowski sum array is empty.

### Input

- `cms` -- cached Minkowski sum

### Output

`true` iff any of the wrapped sets are empty.

### Notes

Forgotten sets cannot be checked anymore.
Normally they should not have been empty because otherwise the support-vector
query would have crashed before.
In that case, the cached Minkowski sum should not be used further.
"""
function isempty(cms::CachedMinkowskiSumArray)
    return any(isempty, array(cms))
end
