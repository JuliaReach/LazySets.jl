"""
    center(B::Ball2)

Return the center of a ball in the 2-norm.

### Input

- `B` -- ball in the 2-norm

### Output

The center of the ball in the 2-norm.
"""
function center(B::Ball2)
    return B.center
end
