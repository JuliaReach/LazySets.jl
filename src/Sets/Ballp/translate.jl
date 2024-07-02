"""
    translate!(B::Ballp, v::AbstractVector)

Translate (i.e., shift) a ball in the p-norm by a given vector, in-place.

### Input

- `B` -- ball in the p-norm
- `v` -- translation vector

### Output

The ball `B` translated by `v`.

### Algorithm

We add the vector to the center of the ball.

### Notes

See also [`translate(::Ballp, ::AbstractVector)`](@ref) for the out-of-place
version.
"""
function translate!(B::Ballp, v::AbstractVector)
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = center(B)
    c .+= v
    return B
end
