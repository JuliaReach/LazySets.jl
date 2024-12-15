"""
    radius_hyperrectangle(x::Interval, i::Int)

Return the box radius of an interval in a given dimension.

### Input

- `x` -- interval
- `i` -- dimension index (must be `1`)

### Output

The box radius in the given dimension.
"""
function radius_hyperrectangle(x::Interval, i::Int)
    @assert i == 1 "an interval has dimension 1, but the index is $i"
    return radius(x)
end

"""
    radius_hyperrectangle(x::Interval)

Return the box radius of an interval in every dimension.

### Input

- `x` -- interval

### Output

The box radius of the interval (a one-dimensional vector).
"""
function radius_hyperrectangle(x::Interval)
    return [radius(x)]
end
