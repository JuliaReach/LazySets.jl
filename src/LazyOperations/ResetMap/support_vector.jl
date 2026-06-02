"""
    σ(d::AbstractVector, rm::ResetMap)

Return a support vector of a reset map.

### Input

- `d`  -- direction
- `rm` -- reset map

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.
"""
@validate function σ(d::AbstractVector, rm::ResetMap)
    N = promote_type(eltype(d), eltype(rm))
    d_reset = copy(d)
    for var in keys(rm.resets)
        d_reset[var] = zero(N)
    end
    return substitute(rm.resets, σ(d_reset, rm.X))
end
