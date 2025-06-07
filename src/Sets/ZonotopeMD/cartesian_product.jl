"""
    cartesian_product(Z1::ZonotopeMD, Z2::ZonotopeMD)

Return the Cartesian product of two structured zonotopes.

### Input

- `Z1` -- structured zonotope
- `Z2` -- structured zonotope

### Output

A new `ZonotopeMD` representing the Cartesian product `Z1 × Z2`.
"""
function cartesian_product(Z1::ZonotopeMD, Z2::ZonotopeMD)
    c = vcat(Z1.center, Z2.center)
    M = blockdiag(sparse(Z1.M), sparse(Z2.M))
    d = vcat(Z1.d, Z2.d)
    return ZonotopeMD(c, M, d)
end
