"""
    AbstractPolynomialZonotope{N} <: LazySet{N}

Abstract type for polynomial zonotope sets.

### Notes

Polynomial zonotopes are in general non-convex. They are always bounded.

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(LazySets.AbstractPolynomialZonotope)
1-element Vector{Any}:
 DensePolynomialZonotope
```
"""
abstract type AbstractPolynomialZonotope{N} end

isconvextype(::Type{<:AbstractPolynomialZonotope}) = false
isboundedtype(::Type{<:AbstractPolynomialZonotope}) = true
