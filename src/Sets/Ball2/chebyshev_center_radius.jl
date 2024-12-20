"""
    chebyshev_center_radius(B::Ball2; [kwargs]...)

Compute a [Chebyshev center](https://en.wikipedia.org/wiki/Chebyshev_center)
and the corresponding radius of a ball in the 2-norm.

### Input

- `B`      -- ball in the 2-norm
- `kwargs` -- further keyword arguments (ignored)

### Output

The pair `(c, r)` where `c` is the Chebyshev center of `B` and `r` is the radius
of the largest Euclidean ball with center `c` enclosed by `B`.

### Notes

The Chebyshev center of a ball in the 2-norm is just the center of the ball.
"""
function chebyshev_center_radius(B::Ball2; kwargs...)
    return B.center, B.radius
end
