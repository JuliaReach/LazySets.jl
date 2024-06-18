"""
    reflect(x::Interval)

Concrete reflection of a interval `x`, resulting in the reflected set `-x`.

### Input

- `x` -- interval

### Output

The `Interval` representing `-x`.

### Algorithm

If ``x = [a, b]``, then ``-x = [-b, -a]``.
"""
function reflect(x::Interval)
    return Interval(-max(x), -min(x))
end
