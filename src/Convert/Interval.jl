"""
    convert(::Type{Interval}, r::Rectification{N, IN}) where {N, IN<:Interval}

Convert a rectification of an interval to an interval.

### Input

- `Interval` -- target type
- `r`        -- rectification of an interval

### Output

An `Interval`.
"""
function convert(::Type{Interval},
                 r::Rectification{N,IN}) where {N,IN<:Interval}
    return Interval(rectify([min(r.X), max(r.X)]))
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
