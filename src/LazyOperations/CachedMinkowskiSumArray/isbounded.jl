"""
    isbounded(cms::CachedMinkowskiSumArray)

Check whether a cached Minkowski sum is bounded.

### Input

- `cms` -- cached Minkowski sum

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(cms::CachedMinkowskiSumArray)
    return all(isbounded, cms.array)
end
