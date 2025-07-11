"""
    norm(MZ::MatrixZonotope, p::Real=Inf)

Compute the operator ``p-norm`` of a matrix zonotope.

### Definition

For a matrix zonotope `\\mathcal{A}``, its ``p``-norm is defined as

```
\\|\\mathcal{A}\\|_p = \\sup_{A \\in \\mathcal{A}} \\|A\\|_p
```

where \``\\|A\\|_p\`` denotes the induced matrix norm.

### Input

- `MZ` -- matrix zonotope set
- `p`  -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.
"""
function norm(MZ::MatrixZonotope, p::Real=Inf)
    return _matrixzonotope_norm(MZ, p, norm)
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
    Gs = generators(MZ)
    n, d = size(C)
    k = ngens(MZ)
    Zmat = Matrix{N}(undef, d, k)
    best = -Inf
    @inbounds for i in 1:n
        c = @view(C[i, :])
        for j in 1:k
            Zmat[:, j] = @view(Gs[j][i, :])
        end
        Z = Zonotope(vec(c), Zmat)
        best = max(best, norm_fn(Z))
    end
    return best
end

function _matrixzonotope_norm(MZ::MatrixZonotope, p::Real, norm_fn::Function)
    dir = p == 1 ? transpose(MZ) : p == Inf ? MZ :
          throw(ArgumentError("the norm for p=$p has not been implemented"))
    return _rowwise_zonotope_norm(dir, norm_fn)
end