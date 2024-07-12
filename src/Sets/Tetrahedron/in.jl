"""
    ∈(x::AbstractVector, T::Tetrahedron)

Check whether a given point is contained in a tetrahedron.

### Input

- `x` -- point/vector
- `P` -- tetrahedron in vertex representation

### Output

`true` iff ``x ∈ P``.

### Algorithm

For each plane of the tetrahedron, we check if the point `x` is on the same side as the remaining vertex.
We need to check this for each plane [1].

[1] https://stackoverflow.com/questions/25179693/how-to-check-whether-the-point-is-in-the-tetrahedron-or-not
"""
function ∈(x::AbstractVector, T::Tetrahedron)
    v = T.vertices
    return same_side(v[1], v[2], v[3], v[4], x) &&
           same_side(v[4], v[1], v[2], v[3], x) &&
           same_side(v[3], v[4], v[1], v[2], x) &&
           same_side(v[2], v[3], v[4], v[1], x)
end

# Return `true` iff point `x` lies in the same half-space as `v4` with respect to the hyperplane `H` determined by `v1`, `v2` and `v3`.
function same_side(v1, v2, v3, v4, x)
    normal = cross(v2 - v1, v3 - v1)
    dotx = dot(normal, x - v1)
    if isapproxzero(dotx)
        return true
    end
    dotv4 = dot(normal, v4 - v1)
    # If the point `x` lies in `H` then `dotx` is zero.
    return signbit(dotv4) == signbit(dotx)
end
