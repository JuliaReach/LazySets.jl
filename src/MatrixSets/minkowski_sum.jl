"""
    minkowski_sum(A::MatrixZonotope, B::MatrixZonotope)

Compute the minkowski sum of two matrix zonotopes.

# Input

- `A` -- a matrix zonotope
- `B` -- a matrix zonotope

# Output

A matrix zonotope representing the minkowski sum between matrix zonotopes.
"""
function minkowski_sum(A::MatrixZonotope, B::MatrixZonotope)
    @assert size(A) == size(B) "cannot sum a matrix zonotope of size $(size(A)) " *
                               "with one of size $(size(B)) "

    c = center(A) + center(B)

    if eltype(A) == eltype(B)
        gensA = generators(A)
        gensB = generators(B)
    else
        gensA = [convert(typeof(c), Ai) for Ai in generators(A)]
        gensB = [convert(typeof(c), Bi) for Bi in generators(B)]
    end

    gens = vcat(gensA, gensB)
    return MatrixZonotope(c, gens)
end
