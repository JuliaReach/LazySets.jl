"""
    ispolytopictype(T::Type{<:LazySet})

Check whether a set type only represents polytopic sets.

### Input

- `T` -- set type

### Output

`true` iff the set type only represents polytopic sets.

### Notes

See also [`ispolytopic(::LazySet)`](@ref).
"""
function ispolytopictype(::Type{<:LazySet}) end
