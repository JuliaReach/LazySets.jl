"""
# Extended help

    σ(d::AbstractVector, U::Universe)

### Output

A vector with infinity values, except in dimensions where the direction is zero.
"""
@validate function σ(d::AbstractVector, U::Universe)
    N = promote_type(eltype(d), eltype(U))
    vec = Vector{N}(undef, length(d))
    @inbounds for (i, vi) in enumerate(d)
        vec[i] = iszero(vi) ? vi : vi > zero(N) ? N(Inf) : N(-Inf)
    end
    return vec
end
