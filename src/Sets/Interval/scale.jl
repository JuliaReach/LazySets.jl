"""
    scale(α::Real, x::Interval)

Concrete scaling of an interval.

### Input

- `α` -- scalar
- `x` -- interval

### Output

The interval obtained by scaling the given interval.
"""
function scale(α::Real, x::Interval)
    return Interval(α * x.dat)
end
