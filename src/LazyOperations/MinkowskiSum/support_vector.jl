"""
    σ(d::AbstractVector, ms::MinkowskiSum)

Return a support vector of a Minkowski sum of two sets.

### Input

- `d`  -- direction
- `ms` -- Minkowski sum of two sets

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.

### Algorithm

A valid support vector in direction ``d`` of the Minkowski sum of two sets ``X``
and ``Y`` is the sum of the support vectors of ``X`` and ``Y`` in direction
``d``.
"""
@validate function σ(d::AbstractVector, ms::MinkowskiSum)
    return σ(d, ms.X) + σ(d, ms.Y)
end
