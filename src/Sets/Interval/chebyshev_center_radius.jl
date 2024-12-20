"""
    chebyshev_center_radius(x::Interval; [kwargs]...)

Compute a [Chebyshev center](https://en.wikipedia.org/wiki/Chebyshev_center)
and the corresponding radius of an interval.

### Input

- `x`      -- interval
- `kwargs` -- further keyword arguments (ignored)

### Output

The pair `(c, r)` where `c` is the Chebyshev center of `x` and `r` is the radius
of the largest Euclidean ball with center `c` enclosed by `x`.

### Notes

The Chebyshev center of an interval is just the center of the interval.
"""
function chebyshev_center_radius(x::Interval; kwargs...)
    return center(x), radius(x)
end
