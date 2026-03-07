"""
    minkowski_sum(A::MatrixZonotope, B::MatrixZonotope)

Compute the Minkowski sum of two matrix zonotopes.

# Input

- `A` -- matrix zonotope
- `B` -- matrix zonotope

# Output

A matrix zonotope representing the Minkowski sum.
"""
function minkowski_sum(A::MatrixZonotope, B::MatrixZonotope)
    @assert size(A) == size(B) "cannot sum matrix zonotopes of size $(size(A)) and $(size(B))"

    c = center(A) + center(B)
    gens = vcat(generators(A), generators(B))
    return MatrixZonotope(c, gens)
end
