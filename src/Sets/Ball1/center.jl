"""
    center(B::Ball1)

Return the center of a ball in the 1-norm.

### Input

- `B` -- ball in the 1-norm

### Output

The center of the ball in the 1-norm.
"""
function center(B::Ball1)
    return B.center
end
