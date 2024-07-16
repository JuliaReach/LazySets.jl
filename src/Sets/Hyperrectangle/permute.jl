"""
    permute(H::Hyperrectangle, p::AbstractVector{Int})

Permute the dimensions according to a permutation vector.

### Input

- `H` -- hyperrectangle
- `p` -- permutation vector

### Output

A permuted hyperrectangle.
"""
function permute(H::Hyperrectangle, p::AbstractVector{Int})
    c = H.center[p]
    r = H.radius[p]
    return Hyperrectangle(c, r)
end
