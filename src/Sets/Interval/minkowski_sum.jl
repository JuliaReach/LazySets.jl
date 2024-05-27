"""
    minkowski_sum(x::Interval, y::Interval)

Concrete Minkowski sum of a pair of intervals.

### Input

- `x` -- interval
- `y` -- interval

### Output

An `Interval` corresponding to the Minkowski sum of `x` and `y`.

### Algorithm

The implementation takes the sum of `x` and `y` following the rules of interval
arithmetic.
"""
function minkowski_sum(x::Interval, y::Interval)
    return Interval(x.dat + y.dat)
end
