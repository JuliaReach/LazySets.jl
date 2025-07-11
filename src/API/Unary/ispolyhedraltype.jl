"""
    ispolyhedraltype(T::Type{<:LazySet})

Check whether a set type only represents polyhedral sets.

### Input

- `T` -- set type

### Output

`true` iff the set type only represents polyhedral sets.

### Notes

See also [`ispolyhedral(::LazySet)`](@ref).
"""
function ispolyhedraltype(::Type{<:LazySet}) end
