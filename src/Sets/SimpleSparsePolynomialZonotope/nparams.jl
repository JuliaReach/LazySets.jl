"""
    nparams(P::SimpleSparsePolynomialZonotope)

Return the number of parameters in the polynomial representation of a simple
sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The number of parameters in the polynomial representation of P.

### Notes

This number corresponds to the number of rows in the exponent matrix ``E`` (`p`
in the mathematical set definition).

### Examples

```jldoctest
julia> S = SimpleSparsePolynomialZonotope([2.0, 0], [1 2;2 2.], [1 4;1 2])
SimpleSparsePolynomialZonotope{Float64, Vector{Float64}, Matrix{Float64}, Matrix{Int64}}([2.0, 0.0], [1.0 2.0; 2.0 2.0], [1 4; 1 2])

julia> nparams(S)
2
```
"""
nparams(P::SSPZ) = size(P.E, 1)
