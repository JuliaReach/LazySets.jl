"""
    diameter(x::Interval, [p]::Real=Inf)

Compute the diameter of an interval, defined as ``‖b - a‖`` in the
``p`-norm, where ``a`` (resp. ``b``) are the minimum (resp. maximum) value of
the interval.

### Input

- `x` -- interval
- `p` -- (optional, default: `Inf`) norm (ignored)

### Output

A real number representing the diameter.

### Notes

In one dimension, all p-norms are identical.
"""
function diameter(x::Interval, ::Real=Inf)
    return max(x) - min(x)
end
