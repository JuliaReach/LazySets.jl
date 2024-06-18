"""
    high(x::Interval)

Return the higher coordinate of an interval set.

### Input

- `x` -- interval

### Output

A vector with the higher coordinate of the interval.
"""
function high(x::Interval)
    return [x.dat.hi]
end
