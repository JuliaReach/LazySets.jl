"""
    permute(U::Universe, p::AbstractVector{Int})

Permute the dimensions according to a permutation vector.

### Input

- `U` -- universe
- `p` -- permutation vector

### Output

The same universe.
"""
function permute(U::Universe, p::AbstractVector{Int})
    return U
end
