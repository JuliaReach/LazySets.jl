"""
    volume(∅::EmptySet{N}) where {N}

Return the volume of an empty set.

### Input

- `∅` -- empty set

### Output

``0``.
"""
function volume(::EmptySet{N}) where {N}
    return zero(N)
end
