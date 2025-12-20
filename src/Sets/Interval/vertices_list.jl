"""
# Extended help

    vertices_list(X::Interval; kwargs...)

### Output

The list of vertices of the interval, which are two one-dimensional vectors,
or just one if the interval is degenerate (the endpoints match within the
working tolerance).
"""
@validate function vertices_list(X::Interval; kwargs...)
    a = min(X)
    b = max(X)
    return _isapprox(a, b) ? [[a]] : [[a], [b]]
end
