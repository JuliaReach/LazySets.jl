"""
    center(x::Interval)

Return the center of an interval.

### Input

- `x` -- interval

### Output

The center, or midpoint, of `x`.
"""
function center(x::Interval)
    return [_center(x)]
end

"""
    center(H::Interval, i::Int)

Return the center along a given dimension of a interval.

### Input

- `x` -- interval
- `i` -- dimension of interest

### Output

The center along a given dimension of the interval.
"""
@inline function center(x::Interval, i::Int)
    @assert i == 1 "an interval has dimension 1, but the index is $i"
    return _center(x)
end

# returns a number, not a vector
function _center(x::Interval)
    return IA.mid(x.dat)
end
