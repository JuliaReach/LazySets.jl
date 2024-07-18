"""
    translate!(Z::Zonotope, v::AbstractVector)

Translate (i.e., shift) a zonotope by a given vector in-place.

### Input

- `Z` -- zonotope
- `v` -- translation vector

### Output

A translated zonotope.

### Notes

See also [`LazySets.API.translate(Z::LazySet, v::AbstractVector)`](@ref) for the
out-of-place version.

### Algorithm

We add the translation vector to the center of the zonotope.
"""
function translate!(Z::Zonotope, v::AbstractVector)
    @assert length(v) == dim(Z) "cannot translate a $(dim(Z))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    Z.center .+= v
    return Z
end
