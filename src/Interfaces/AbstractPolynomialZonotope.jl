"""
    AbstractPolynomialZonotope{N} <: LazySet{N}

Abstract type for polynomial zonotope sets.

### Notes

Polynomial zonotopes are in general non-convex. They are always bounded.

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(LazySets.AbstractPolynomialZonotope)
2-element Vector{Any}:
 DensePolynomialZonotope
 SimpleSparsePolynomialZonotope
```
"""
abstract type AbstractPolynomialZonotope{N} end

isconvextype(::Type{<:AbstractPolynomialZonotope}) = false
isboundedtype(::Type{<:AbstractPolynomialZonotope}) = true
