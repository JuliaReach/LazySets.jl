"""
    translate!(B::Ball2, v::AbstractVector)

Translate (i.e., shift) a ball in the 2-norm by the given vector, in-place.

### Input

- `B` -- ball in the 2-norm
- `v` -- translation vector

### Output

The ball `B` translated by `v`.

### Algorithm

We add the vector to the center of the ball.

### Notes

See also [`translate(::Ball2, ::AbstractVector)`](@ref) for the out-of-place version.
"""
function translate!(B::Ball2, v::AbstractVector)
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = B.center
    c .+= v
    return B
end
