
"""
    cartesian_product(Z1::ZonotopeMD, Z2::ZonotopeMD)

Return the Cartesian product of two zonotopes in normal form (`ZonotopeMD`).

### Input

- `Z1`, `Z2` -- zonotopes in normal form (`ZonotopeMD`)

### Output

A new `ZonotopeMD` representing the Cartesian product `Z1 × Z2`.
"""
function cartesian_product(Z1::ZonotopeMD, Z2::ZonotopeMD)
    c = vcat(Z1.center, Z2.center)
    d = vcat(Z1.d, Z2.d)
    M = blockdiag(sparse(Z1.M), sparse(Z2.M))
    return ZonotopeMD(c, M, d)
end
