"""
    ∈(x::AbstractVector, Z::ZeroSet)

Check whether a given point is contained in a zero set.

### Input

- `x` -- point/vector
- `Z` -- zero set

### Output

`true` iff ``x ∈ Z``.

### Examples

```jldoctest
julia> Z = ZeroSet(2);

julia> [1.0, 0.0] ∈ Z
false
julia> [0.0, 0.0] ∈ Z
true
```
"""
function ∈(x::AbstractVector, Z::ZeroSet)
    @assert length(x) == dim(Z) "a $(length(x))-dimensional vector is " *
                                "incompatible with a $(dim(Z))-dimensional set"
    return iszero(x)
end
