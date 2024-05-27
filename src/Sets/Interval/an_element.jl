"""
    an_element(x::Interval)

Return some element of an interval.

### Input

- `x` -- interval

### Output

The left border (`low(x)`) of the interval.
"""
function an_element(x::Interval)
    return low(x)
end
