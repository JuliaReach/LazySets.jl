"""
    low(Z::Zonotope, i::Int)

Return the lower coordinate of a zonotope in a given dimension.

### Input

- `Z` -- zonotope
- `i` -- dimension of interest

### Output

The lower coordinate of the zonotope in the given dimension.
"""
function low(Z::Zonotope, i::Int)
    G = genmat(Z)
    v = center(Z, i)
    @inbounds for j in 1:ngens(Z)
        v -= abs(G[i, j])
    end
    return v
end
