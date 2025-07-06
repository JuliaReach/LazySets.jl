"""
    ngens(X::Interval)

Return the number of generators of an interval.

### Input

- `X` -- interval

### Output

The number of generators.

### Algorithm

An interval has either one generator, or zero generators if it is a degenerated
interval of diameter zero.
"""
function ngens(X::Interval)
    return _isapprox(min(X), max(X)) ? 0 : 1
end
