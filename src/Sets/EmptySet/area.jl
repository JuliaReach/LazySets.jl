"""
    area(∅::EmptySet)

Return the area of an empty set.

### Input

- `∅` -- empty set

### Output

``0``.
"""
function area(∅::EmptySet)
    N = eltype(∅)
    return zero(N)
end
