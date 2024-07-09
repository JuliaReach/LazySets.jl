"""
    translate!(E::Ellipsoid, v::AbstractVector)

Translate (i.e., shift) an ellipsoid by a given vector, in-place.

### Input

- `E` -- ellipsoid
- `v` -- translation vector

### Output

The ellipsoid `E` translated by `v`.

### Notes

See also [`translate(::Ellipsoid, ::AbstractVector)`](@ref) for the out-of-place
version.

### Algorithm

We add the vector to the center of the ellipsoid.
"""
function translate!(E::Ellipsoid, v::AbstractVector)
    @assert length(v) == dim(E) "cannot translate a $(dim(E))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = E.center
    c .+= v
    return E
end
