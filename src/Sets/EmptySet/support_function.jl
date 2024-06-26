"""
    ρ(d::AbstractVector, ∅::EmptySet)

Evaluate the support function of an empty set in a given direction.

### Input

- `d` -- direction
- `∅` -- empty set

### Output

An error.
"""
function ρ(::AbstractVector, ::EmptySet)
    return error("the support function of an empty set is undefined")
end
