"""
    high(Z::Zonotope, i::Int)

Return the higher coordinate of a zonotope in a given dimension.

### Input

- `Z` -- zonotope
- `i` -- dimension of interest

### Output

The higher coordinate of the zonotope in the given dimension.
"""
function high(Z::Zonotope, i::Int)
    G = genmat(Z)
    v = center(Z, i)
    @inbounds for j in 1:ngens(Z)
        v += abs(G[i, j])
    end
    return v
end
