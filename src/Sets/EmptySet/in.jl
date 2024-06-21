"""
    ∈(x::AbstractVector, ∅::EmptySet)

Check whether a given point is contained in an empty set.

### Input

- `x` -- point/vector
- `∅` -- empty set

### Output

`false`.

### Examples

```jldoctest
julia> [1.0, 0.0] ∈ ∅(2)
false
```
"""
function ∈(::AbstractVector, ::EmptySet)
    return false
end
