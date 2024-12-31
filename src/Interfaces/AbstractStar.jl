"""
    AbstractStar{N} <: LazySet{N}

Abstract supertype for all star set types.

### Notes

A set ``X`` is star-like (also known as generalized star) if it can be
represented by a center ``x₀ ∈ ℝ^n`` and ``m`` vectors ``v₁, …, vₘ``
forming the basis, and a predicate ``P : ℝ^n → \\{⊤, ⊥\\}`` such that

```math
    X = \\{x ∈ ℝ^n : x = x₀ + ∑_{i=1}^m α_i v_i,~~\\textrm{s.t. } P(α) = ⊤ \\}.
```
"""
abstract type AbstractStar{N} <: LazySet{N} end

isoperationtype(::Type{<:AbstractStar}) = false
