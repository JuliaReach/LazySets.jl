"""
    isbounded(X::LazySet)

Check whether a set is bounded.

### Input

- `X` -- set

### Output

`true` iff the set is bounded.

### Notes

See also [`isboundedtype(::Type{<:LazySet})`](@ref).
"""
function isbounded(::LazySet) end
