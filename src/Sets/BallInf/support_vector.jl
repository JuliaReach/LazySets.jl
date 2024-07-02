"""
    σ(d::AbstractVector, B::BallInf)

Return the support vector of a ball in the infinity norm in the given direction.

### Input

- `d` -- direction
- `B` -- ball in the infinity norm

### Output

The support vector in the given direction.
If the direction has norm zero, the center of the ball is returned.
"""
function σ(d::AbstractVector, B::BallInf)
    @assert length(d) == dim(B) "a $(length(d))-dimensional vector is " *
                                "incompatible with a $(dim(B))-dimensional set"
    return center(B) .+ sign.(d) .* B.radius
end

# special case for SingleEntryVector
function σ(d::SingleEntryVector, B::BallInf)
    return _σ_sev_hyperrectangle(d, B)
end
