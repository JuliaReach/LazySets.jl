"""
    norm(MZ::MatrixZonotope, p::Real=Inf)

Compute the operator ``p``-norm of a matrix zonotope.

### Definition

For a matrix zonotope `\\mathcal{A}``, its ``p``-norm is defined as

```
\\|\\mathcal{A}\\|_p = \\sup_{A \\in \\mathcal{A}} \\|A\\|_p
```

where ``\\|A\\|_p`` denotes the induced matrix norm.

### Input

- `MZ` -- matrix zonotope set
- `p`  -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.
"""
function norm(MZ::MatrixZonotope, p::Real=Inf)
    if p == 1
        return _rowwise_zonotope_norm(transpose(MZ), norm)
    elseif p == Inf
        return _rowwise_zonotope_norm(MZ, norm)
    else
        throw(ArgumentError("the norm for p=$p has not been implemented"))
    end
end

"""
    _rowwise_zonotope_norm(MZ::MatrixZonotope{N}, norm_fn::Function) where {N}

Compute the induced matrix norm of a matrix zonotope by reducing to row-wise zonotope ``ℓ₁`` norms.

### Input

- `MZ` -- matrix zonotope set
- `norm_fn` -- a function to approximate the zonotope ``ℓ₁`` norm

### Output

The induced matrix ``p``-norm of the matrix zonotope.

### Algorithm

For each row index `i = 1, ..., n`, we construct a zonotope with center given 
by the `i`-th row of the center matrix and as generators the `i`-th row of each generator matrix.
The norm of this zonotope is then computed using the provided `norm_fn`.
The final result is the maximum of these `n` row-wise zonotope norms.
"""
function _rowwise_zonotope_norm(MZ::MatrixZonotope{N}, norm_fn::Function) where {N}
    C = center(MZ)
    n, d = size(C)
    k = ngens(MZ)
    Gmat = Matrix{N}(undef, d, k)
    best = -Inf
    @inbounds for i in 1:n
        c = @view(C[i, :])
        for (j, G) in enumerate(generators(MZ))
            Gmat[:, j] = @view(G[i, :])
        end
        Z = Zonotope(c, Gmat)
        best = max(best, norm_fn(Z, 1))
    end
    return best
end
