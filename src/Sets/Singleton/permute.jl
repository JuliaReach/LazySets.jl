"""
    permute(S::Singleton, p::AbstractVector{Int})

Permute the dimensions according to a permutation vector.

### Input

- `S` -- singleton
- `p` -- permutation vector

### Output

A new `Singleton` with the permuted dimensions.
"""
function permute(S::Singleton, p::AbstractVector{Int})
    e = S.element[p]
    return Singleton(e)
end
