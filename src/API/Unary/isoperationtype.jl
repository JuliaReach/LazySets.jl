"""
    isoperationtype(T::Type{<:LazySet})

Check whether a set type only represents (lazy) set operations.

### Input

- `T` -- set type

### Output

`true` iff the set type only represents set operations.

### Notes

See also [`isoperation(::LazySet)`](@ref).
"""
function isoperationtype(::Type{<:LazySet}) end
