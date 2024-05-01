"""
    concretize(X::LazySet)

Construct a concrete representation of a (possibly lazy) set.

### Input

- `X` -- set

### Output

A concrete representation of `X` (as far as possible).

### Notes

Since not every lazy set has a concrete set representation in this library, the
result may still be partially lazy.
"""
function concretize(::LazySet) end
