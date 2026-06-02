"""
    ρ(d::AbstractVector, rm::ResetMap)

Evaluate the support function of a reset map.

### Input

- `d`  -- direction
- `rm` -- reset map

### Output

The evaluation of the support function in the given direction.

### Notes

We use the usual dot-product definition, but for unbounded sets we redefine the
product between ``0`` and ``±∞`` as ``0``; Julia returns `NaN` here.

```jldoctest
julia> Inf * 0.0
NaN
```

See the discussion [here](https://math.stackexchange.com/q/28940).
"""
@validate function ρ(d::AbstractVector, rm::ResetMap)
    return dot_zero(d, σ(d, rm))
end
