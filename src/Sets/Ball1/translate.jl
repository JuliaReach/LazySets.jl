"""
    translate!(B::Ball1, v::AbstractVector)

Translate (i.e., shift) a ball in the 1-norm by the given vector, in-place.

### Input

- `B` -- ball in the 1-norm
- `v` -- translation vector

### Output

The in-place translated ball in the 1-norm.

### Algorithm

We add the vector to the center of the ball.

### Notes

See also [`translate(::Ball1, ::AbstractVector)`](@ref) for the out-of-place version.
"""
function translate!(B::Ball1, v::AbstractVector)
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = B.center
    c .+= v
    return B
end
