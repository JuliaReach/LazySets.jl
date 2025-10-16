"""
    ispolytopic(X::LazySet)

Check whether a set is polytopic.

### Input

- `X` -- set

### Output

`true` only if the set is polytopic.

### Notes

The answer is conservative, i.e., may sometimes be `false` even if the set is
polytopic.

A set is polytopic if it is both polyhedral and bounded.
"""
function ispolytopic(::LazySet) end
