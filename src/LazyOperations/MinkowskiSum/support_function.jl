"""
    ρ(d::AbstractVector, ms::MinkowskiSum)

Evaluate the support function of a Minkowski sum of two sets.

### Input

- `d`  -- direction
- `ms` -- Minkowski sum of two sets

### Output

The evaluation of the support function in the given direction.

### Algorithm

The support function in direction ``d`` of the Minkowski sum of two sets ``X``
and ``Y`` is the sum of the support functions of ``X`` and ``Y`` in direction
``d``.
"""
@validate function ρ(d::AbstractVector, ms::MinkowskiSum)
    return ρ(d, ms.X) + ρ(d, ms.Y)
end
