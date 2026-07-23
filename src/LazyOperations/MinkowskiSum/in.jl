"""
    in(x::AbstractVector, ms::MinkowskiSum{N,<:AbstractSingleton}) where {N}

Check whether a given point is contained in the Minkowski sum of a singleton
and another set.

### Input

- `x`  -- point/vector
- `ms` -- Minkowski sum of a singleton and another set

### Output

`true` iff ``x ∈ ms``.

### Algorithm

Note that ``x ∈ (S ⊕ P)``, where ``S = \\{s\\}``  is a singleton set and
``P`` is a set, if and only if ``(x-s) ∈ P``.
"""
@validate function in(x::AbstractVector, ms::MinkowskiSum{N,<:AbstractSingleton}) where {N}
    return _in_singleton_msum(x, ms.X, ms.Y)
end

# symmetric method
@validate function in(x::AbstractVector,
                      ms::MinkowskiSum{N,<:LazySet,<:AbstractSingleton}) where {N}
    return _in_singleton_msum(x, ms.Y, ms.X)
end

# disambiguation
@validate function in(x::AbstractVector,
                      ms::MinkowskiSum{N,<:AbstractSingleton,<:AbstractSingleton}) where {N}
    return _in_singleton_msum(x, ms.X, ms.Y)
end

@inline _in_singleton_msum(x, X, Y) = (x - center(X)) ∈ Y
