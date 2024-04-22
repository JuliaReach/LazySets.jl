"""
    isoperation(X::LazySet)

Check whether a set is an instance of a (lazy) set operation.

### Input

- `X` -- set

### Output

`true` iff `X` is an instance of a set-based operation.

### Notes

See also [`isoperationtype(::Type{<:LazySet})`](@ref).
"""
function isoperation(::LazySet) end
