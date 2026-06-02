"""
    σ(d::AbstractVector, lm::LinearMap)

Return a support vector of the linear map.

### Input

- `d`  -- direction
- `lm` -- linear map

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``L = M⋅S``, where ``M`` is a matrix and ``S`` is a set, it follows that
``σ(d, L) = M⋅σ(M^T d, S)`` for any direction ``d``.
"""
@validate function σ(d::AbstractVector, lm::LinearMap)
    return _σ_linear_map(d, lm.M, lm.X)
end

function _σ_linear_map(d::AbstractVector, M::AbstractMatrix, X::LazySet)
    return M * σ(At_mul_B(M, d), X)
end
