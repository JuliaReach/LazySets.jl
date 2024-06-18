"""
    low(x::Interval)

Return the lower coordinate of an interval set.

### Input

- `x` -- interval

### Output

A vector with the lower coordinate of the interval.
"""
function low(x::Interval)
    return [x.dat.lo]
end
