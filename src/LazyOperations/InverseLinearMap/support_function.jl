"""
    ρ(d::AbstractVector, ilm::InverseLinearMap)

Evaluate the support function of the inverse linear map.

### Input

- `d`      -- direction
- `ilm`    -- inverse linear map

### Output

The evaluation of the support function in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``L = M^{-1}⋅X``, where ``M`` is a matrix and ``X`` is a set, it follows
that ``ρ(d, L) = ρ((M^T)^{-1} d, X)`` for any direction ``d``.
"""
@validate function ρ(d::AbstractVector, ilm::InverseLinearMap)
    y = transpose(ilm.M) \ d
    return ρ(y, ilm.X)
end
