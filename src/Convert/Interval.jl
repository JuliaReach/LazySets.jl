"""
    convert(::Type{Interval}, R::Rectification{N,<:Interval}) where {N}

Convert a rectification of an interval to an interval.

### Input

- `Interval` -- target type
- `R`        -- rectification of an interval

### Output

An `Interval`.
"""
function convert(::Type{Interval}, R::Rectification{N,<:Interval}) where {N}
    return rectify(R.X)
end

"""
    convert(::Type{Interval}, ms::MinkowskiSum{N, IT, IT}) where {N, IT<:Interval}

Convert the Minkowski sum of two intervals to an interval.

### Input

- `Interval` -- target type
- `ms`       -- Minkowski sum of two intervals

### Output

An interval.
"""
function convert(::Type{Interval}, ms::MinkowskiSum{N,IT,IT}) where {N,IT<:Interval}
    return concretize(ms)
end
