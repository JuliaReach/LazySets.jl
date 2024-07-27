"""
    ispolyhedral(X::LazySet)

Check whether a set is polyhedral.

### Input

- `X` -- set

### Output

`true` only if the set is polyhedral.

### Notes

The answer is conservative, i.e., may sometimes be `false` even if the set is
polyhedral.
"""
function ispolyhedral(::LazySet) end
