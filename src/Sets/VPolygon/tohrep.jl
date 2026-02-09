"""
    tohrep(P::VPolygon, ::Type{HPOLYGON}=HPolygon) where {HPOLYGON<:AbstractHPolygon}

Build a constraint representation of the given polygon.

### Input

- `P`        -- polygon in vertex representation
- `HPOLYGON` -- (optional, default: `HPolygon`) type of target polygon

### Output

An `HPOLYGON` (polygon in constraint representation).
"""
function tohrep(P::VPolygon, ::Type{HPOLYGON}=HPolygon) where {HPOLYGON<:AbstractHPolygon}
    clist = constraints_list(P)
    return HPOLYGON(clist)
end
