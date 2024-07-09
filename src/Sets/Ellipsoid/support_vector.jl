"""
    σ(d::AbstractVector, E::Ellipsoid)

Return the support vector of an ellipsoid in a given direction.

### Input

- `d` -- direction
- `E` -- ellipsoid

### Output

The support vector in the given direction.

### Algorithm

Let ``E`` be an ellipsoid of center ``c`` and shape matrix ``Q = BB^\\mathrm{T}``.
Its support vector along direction ``d`` can be deduced from that of the unit
Euclidean ball ``\\mathcal{B}_2`` using the algebraic relations for the support
vector,

```math
σ_{B\\mathcal{B}_2 ⊕ c}(d) = c + Bσ_{\\mathcal{B}_2} (B^\\mathrm{T} d)
= c + \\dfrac{Qd}{\\sqrt{d^\\mathrm{T}Q d}}.
```
"""
function σ(d::AbstractVector, E::Ellipsoid)
    if iszero(norm(d, 2))
        return E.center
    end
    Qd = E.shape_matrix * d
    return E.center .+ Qd ./ sqrt(dot(d, Qd))
end
