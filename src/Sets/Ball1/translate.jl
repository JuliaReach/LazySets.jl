"""
# Extended help

    translate!(B::Ball1, v::AbstractVector)

### Algorithm

We add the vector to the center of the ball.
"""
function translate!(B::Ball1, v::AbstractVector)
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = B.center
    c .+= v
    return B
end
