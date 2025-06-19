"""
    radius_hyperrectangle(X::Interval, i::Int)

Return the box radius of an interval in a given dimension.

### Input

- `X` -- interval
- `i` -- dimension index (must be `1`)

### Output

The box radius in the given dimension.
"""
function radius_hyperrectangle(X::Interval, i::Int)
    @assert i == 1 "an interval has dimension 1, but the index is $i"
    return radius(X)
end

"""
    radius_hyperrectangle(X::Interval)

Return the box radius of an interval in every dimension.

### Input

- `X` -- interval

### Output

The box radius of the interval (a one-dimensional vector).
"""
function radius_hyperrectangle(X::Interval)
    return [radius(X)]
end
