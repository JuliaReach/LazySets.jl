"""
    vertices_list(x::Interval; kwargs...)

Return the list of vertices of an interval.

### Input

- `x` -- interval

### Output

The list of vertices of the interval, which are two one-dimensional vectors,
or just one if the interval is degenerate (the endpoints match within the
working tolerance).
"""
function vertices_list(x::Interval; kwargs...)
    a = min(x)
    b = max(x)
    return _isapprox(a, b) ? [[a]] : [[a], [b]]
end
