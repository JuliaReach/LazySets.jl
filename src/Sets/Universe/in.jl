"""
    ∈(x::AbstractVector, U::Universe)

Check whether a given point is contained in a universe.

### Input

- `x` -- point/vector
- `U` -- universe

### Output

`true`.

### Examples

```jldoctest
julia> [1.0, 0.0] ∈ Universe(2)
true
```
"""
function ∈(x::AbstractVector, U::Universe)
    @assert length(x) == dim(U) "a $(length(x))-dimensional vector is " *
                                "incompatible with a $(dim(U))-dimensional set"
    return true
end
