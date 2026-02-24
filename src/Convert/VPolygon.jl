"""
    convert(::Type{VPolygon}, P::AbstractHPolygon)

Convert a polygon in constraint representation to a polygon in vertex
representation.

### Input

- `VPolygon` -- target type
- `P`        -- polygon in constraint representation

### Output

A polygon in vertex representation.
"""
function convert(::Type{VPolygon}, P::AbstractHPolygon)
    return tovrep(P)
end
