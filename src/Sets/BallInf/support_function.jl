const BALLINF_THRESHOLD_ρ = 30  # threshold value in `ρ`

"""
# Extended help

    ρ(d::AbstractVector, B::BallInf)

### Algorithm

Let ``B`` be a ball in the infinity norm with center ``c`` and radius ``r`` and
let ``d`` be the direction of interest.
For balls with dimensions less than ``$BALLINF_THRESHOLD_ρ``, we use the
implementation for `AbstractHyperrectangle`, tailored to a `BallInf`, which
computes

```math
    ∑_{i=1}^n d_i (c_i + \\textrm{sgn}(d_i) · r)
```

where ``\\textrm{sgn}(α) = 1`` if ``α ≥ 0`` and ``\\textrm{sgn}(α) = -1`` if ``α < 0``.

For balls of higher dimension, we instead exploit that for a support vector
``v = σ(d, B) = c + \\textrm{sgn}(d) · (r, …, r)ᵀ`` we have

```math
    ρ(d, B) = ⟨d, v⟩ = ⟨d, c⟩ + ⟨d, \\textrm{sgn}(d) · (r, …, r)ᵀ⟩ = ⟨d, c⟩ + r · ∑_{i=1}^n |d_i|
```

where ``⟨·, ·⟩`` denotes the dot product.
"""
function ρ(d::AbstractVector, B::BallInf)
    @assert length(d) == dim(B) "a $(length(d))-dimensional vector is " *
                                "incompatible with a $(dim(B))-dimensional set"
    c = center(B)
    r = B.radius
    if length(d) > BALLINF_THRESHOLD_ρ
        # more efficient for higher dimensions
        return dot(d, c) + r * sum(abs, d)
    end
    N = promote_type(eltype(d), eltype(B))
    res = zero(N)
    @inbounds for (i, di) in enumerate(d)
        if di < zero(N)
            res += di * (c[i] - r)
        elseif di > zero(N)
            res += di * (c[i] + r)
        end
    end
    return res
end

# special case for SingleEntryVector
function ρ(d::SingleEntryVector, B::BallInf)
    return _ρ_sev_hyperrectangle(d, B)
end
