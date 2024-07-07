"""
    ρ(d::AbstractVector, Z::ZeroSet)

Evaluate the support function of a zero set in a given direction.

### Input

- `d` -- direction
- `Z` -- zero set

### Output

`0`.
"""
function ρ(d::AbstractVector, Z::ZeroSet)
    @assert length(d) == dim(Z) "a $(length(d))-dimensional vector is " *
                                "incompatible with a $(dim(Z))-dimensional set"
    N = promote_type(eltype(d), eltype(Z))
    return zero(N)
end
