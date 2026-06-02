"""
    ρ(d::AbstractVector, eprojmap::ExponentialProjectionMap;
      [backend]=get_exponential_backend())

Evaluate the support function of a projection of an exponential map.

### Input

- `d`        -- direction
- `eprojmap` -- projection of an exponential map
- `backend`  -- (optional; default: `get_exponential_backend()`) exponentiation
                backend

### Output

Evaluation of the support function in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``S = (L⋅M⋅R)⋅X``, where ``L`` and ``R`` are matrices, ``M`` is a matrix
exponential, and ``X`` is a set, it follows that ``ρ(d, S) = ρ(R^T⋅M^T⋅L^T⋅d, X)``
for any direction ``d``.
"""
@validate function ρ(d::AbstractVector, eprojmap::ExponentialProjectionMap;
                     backend=get_exponential_backend())
    N = promote_type(eltype(d), eltype(eprojmap))
    Lᵀd = transpose(eprojmap.projspmexp.L) * d
    eᴹLᵀd = _expmv(backend, one(N), transpose(eprojmap.projspmexp.spmexp.M), Lᵀd)
    RᵀeᴹLᵀd = At_mul_B(eprojmap.projspmexp.R, eᴹLᵀd)
    return ρ(RᵀeᴹLᵀd, eprojmap.X)
end
