"""
    eltype(T::Type{<:LazySet})

Determine the numeric type of a set type.

### Input

- `T` -- set type

### Output

The numeric type of `T`.
"""
function eltype(::Type{<:LazySet}) end

"""
    eltype(X::LazySet)

Determine the numeric type of a set.

### Input

- `X` -- set

### Output

The numeric type of `X`.
"""
function eltype(::LazySet) end
