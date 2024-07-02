"""
    isflat(B::BallInf)

Determine whether a ball in the infinity norm is flat, i.e., whether its radius
is zero.

### Input

- `B` -- ball in the infinity norm

### Output

`true` iff the ball is flat.

### Notes

For robustness with respect to floating-point inputs, this function relies on
the result of `isapproxzero` applied to the radius of the ball.
Hence, this function depends on the absolute zero tolerance `ABSZTOL`.
"""
function isflat(B::BallInf)
    return isapproxzero(B.radius)
end
