"""
    σ(d::AbstractVector, ilm::InverseLinearMap)

Return a support vector of a inverse linear map.

### Input

- `d`   -- direction
- `ilm` -- inverse linear map

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``L = M^{-1}⋅X``, where ``M`` is a matrix and ``X`` is a set, since
(M^T)^{-1}=(M^{-1})^T, it follows that ``σ(d, L) = M^{-1}⋅σ((M^T)^{-1} d, X)``
for any direction ``d``.
"""
@validate function σ(d::AbstractVector, ilm::InverseLinearMap)
    y = transpose(ilm.M) \ d
    return ilm.M \ σ(y, ilm.X)
end
