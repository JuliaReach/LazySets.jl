"""
    ngens(x::Interval)

Return the number of generators of an interval.

### Input

- `x` -- interval

### Output

The number of generators.

### Algorithm

An interval has either one generator, or zero generators if it is a degenerated
interval of diameter zero.
"""
function ngens(x::Interval)
    return _isapprox(min(x), max(x)) ? 0 : 1
end
