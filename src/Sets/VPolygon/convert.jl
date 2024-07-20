"""
    convert(::Type{VPolygon}, X::LazySet)

Convert a two-dimensional polytopic set to a polygon in vertex representation.

### Input

- `VPolygon` -- target type
- `X`        -- two-dimensional polytopic set

### Output

The 2D set represented as a polygon.

### Algorithm

This method uses `vertices_list`.
"""
function convert(::Type{VPolygon}, X::LazySet)
    @assert dim(X) == 2 "set must be two-dimensional for conversion"
    return VPolygon(vertices_list(X))
end
