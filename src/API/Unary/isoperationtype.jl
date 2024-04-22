"""
    isoperationtype(T::Type{<:LazySet})

Check whether a set type is a (lazy) set operation.

### Input

- `T` -- set type

### Output

`true` iff the set type represents a set operation.

### Notes

See also [`isoperation(::LazySet)`](@ref).
"""
function isoperationtype(::Type{<:LazySet}) end
