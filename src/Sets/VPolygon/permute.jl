"""
    permute(V::VPolygon, p::AbstractVector{Int})

Permute the dimensions according to a permutation vector.

### Input

- `P` -- polygon in vertex representation
- `p` -- permutation vector

### Output

The permuted polygon in vertex representation.
"""
function permute(V::VPolygon, p::AbstractVector{Int})
    return VPolygon([v[p] for v in V.vertices]; apply_convex_hull=true)
end
