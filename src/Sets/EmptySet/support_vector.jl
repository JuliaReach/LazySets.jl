"""
    σ(d::AbstractVector, ∅::EmptySet)

Return the support vector of an empty set.

### Input

- `d` -- direction
- `∅` -- empty set

### Output

An error.
"""
function σ(::AbstractVector, ::EmptySet)
    return error("the support vector of an empty set is undefined")
end
