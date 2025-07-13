"""
# Extended help

    translate!(B::Ball2, v::AbstractVector)

Translate (i.e., shift) a ball in the 2-norm by the given vector, in-place.

### Algorithm

We add the vector to the center of the ball.
"""
@validate function translate!(B::Ball2, v::AbstractVector)
    c = B.center
    c .+= v
    return B
end
