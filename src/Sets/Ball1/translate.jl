"""
# Extended help

    translate!(B::Ball1, v::AbstractVector)

### Algorithm

We add the vector to the center of the ball.
"""
@validate function translate!(B::Ball1, v::AbstractVector)
    c = B.center
    c .+= v
    return B
end
