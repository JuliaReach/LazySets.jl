"""
    norm(MZ::MatrixZonotope, p::Real=Inf)

Compute the operator ``p``-norm of a matrix zonotope.

### Input

- `MZ` -- matrix zonotope
- `p`  -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.

### Notes

For a matrix zonotope `\\mathcal{A}``, its ``p``-norm is defined as

```math
‖\\mathcal{A}‖_p = \\sup_{A \\in \\mathcal{A}} ‖A‖_p
```

where ``‖A‖_p`` denotes the induced matrix norm.
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

- `MZ`      -- matrix zonotope
- `norm_fn` -- function to approximate the zonotope ``ℓ₁`` norm

### Output

The induced matrix ``p``-norm of the matrix zonotope.

### Algorithm

For each row index ``i = 1, ..., n``, we construct a zonotope with center given
by the ``i``-th row of the center matrix and as generators the ``i``-th row of
each generator matrix. The norm of this zonotope is then computed using the provided
`norm_fn`. The final result is the maximum of these ``n`` row-wise zonotope norms.
"""
function _rowwise_zonotope_norm(MZ::MatrixZonotope{N}, norm_fn::Function) where {N}
    A0 = center(MZ)
    n, d = size(A0)
    k = ngens(MZ)
    G = Matrix{N}(undef, d, k)
    best = N(-Inf)
    @inbounds for i in 1:n
        c = @view(A0[i, :])
        for (j, Ai) in enumerate(generators(MZ))
            G[:, j] = @view(Ai[i, :])
        end
        Z = Zonotope(c, G)
        best = max(best, norm_fn(Z, 1))
    end
    return best
end
