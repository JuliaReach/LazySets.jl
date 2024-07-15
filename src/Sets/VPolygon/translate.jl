"""
    translate(P::VPolygon, v::AbstractVector)

Translate (i.e., shift) a polygon in vertex representation by a given vector.

### Input

- `P` -- polygon in vertex representation
- `v` -- translation vector

### Output

A translated polygon in vertex representation.

### Notes

See [`translate!(::VPolygon, ::AbstractVector)`](@ref) for the in-place version.
"""
function translate(P::VPolygon, v::AbstractVector)
    return translate!(deepcopy(P), v)
end

"""
    translate!(P::VPolygon, v::AbstractVector)

Translate (i.e., shift) a polygon in vertex representation by a given vector,
in-place.

### Input

- `P` -- polygon in vertex representation
- `v` -- translation vector

### Output

The polygon translated by the vector.

### Algorithm

We add the vector to each vertex of the polygon.

### Notes

See [`translate(::VPolygon, ::AbstractVector)`](@ref) for the out-of-place
version.
"""
function translate!(P::VPolygon, v::AbstractVector)
    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    for x in P.vertices
        x .+= v
    end
    return P
end
