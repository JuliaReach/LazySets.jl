"""
# Extended help

    translate!(B::Ballp, v::AbstractVector)

### Algorithm

We add the vector to the center of the ball.
"""
@validate function translate!(B::Ballp, v::AbstractVector)
    c = center(B)
    c .+= v
    return B
end
