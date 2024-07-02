"""
    center(B::BallInf)

Return the center of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm

### Output

The center of the ball in the infinity norm.
"""
function center(B::BallInf)
    return B.center
end
