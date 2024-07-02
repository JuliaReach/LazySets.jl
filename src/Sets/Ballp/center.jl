"""
    center(B::Ballp)

Return the center of a ball in the p-norm.

### Input

- `B` -- ball in the p-norm

### Output

The center of the ball in the p-norm.
"""
function center(B::Ballp)
    return B.center
end
