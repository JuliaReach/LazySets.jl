@doc """
    eltype(T::Type{<:LazySet})

Determine the numeric type of a set type.

### Input

- `T` -- set type

### Output

The numeric type of `T`.
"""
eltype(::Type{<:LazySet})
