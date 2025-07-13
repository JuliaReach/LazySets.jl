"""
# Extended help

    translate!(B::BallInf, v::AbstractVector)

### Algorithm

We add the vector to the center of the ball.
"""
@validate function translate!(B::BallInf, v::AbstractVector)
    c = B.center
    c .+= v
    return B
end
