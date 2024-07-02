"""
    translate!(B::BallInf, v::AbstractVector)

Translate (i.e., shift) a ball in the infinity norm by a given vector, in-place.

### Input

- `B` -- ball in the infinity norm
- `v` -- translation vector

### Output

The ball `B` translated by `v`.

### Algorithm

We add the vector to the center of the ball.

### Notes

See also [`translate(::BallInf, ::AbstractVector)`](@ref) for the out-of-place
version.
"""
function translate!(B::BallInf, v::AbstractVector)
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = B.center
    c .+= v
    return B
end
