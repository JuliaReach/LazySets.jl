"""
    σ(d::AbstractVector, B::Ball2)

Return the support vector of a 2-norm ball in the given direction.

### Input

- `d` -- direction
- `B` -- ball in the 2-norm

### Output

The support vector in the given direction.
If the direction has norm zero, the center is returned.

### Notes

Let ``c`` and ``r`` be the center and radius of a ball ``B`` in the 2-norm,
respectively.
For nonzero direction ``d`` we have

```math
σ(d, B) = c + r \\frac{d}{‖d‖_2}.
```

This function requires computing the 2-norm of the input direction, which is
performed in the given precision of the numeric datatype of both the direction
and the set.
Exact inputs are not supported.
"""
function σ(d::AbstractVector, B::Ball2)
    dnorm = norm(d, 2)
    if isapproxzero(dnorm)
        return B.center
    else
        return @. B.center + d * (B.radius / dnorm)
    end
end
