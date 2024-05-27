"""
    σ(d::AbstractVector, x::Interval)

Return the support vector of an interval in a given direction.

### Input

- `d` -- direction
- `x` -- interval

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector, x::Interval)
    @assert length(d) == dim(x) "a $(length(d))-dimensional vector is " *
                                "incompatible with an $(dim(x))-dimensional set"
    N = promote_type(eltype(d), eltype(x))
    return @inbounds d[1] > zero(N) ? high(x) : low(x)
end
