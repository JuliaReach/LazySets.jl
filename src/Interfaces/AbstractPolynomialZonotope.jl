export AbstractPolynomialZonotope

"""
    AbstractPolynomialZonotope{N} <: LazySet{N}

Abstract type for polynomial zonotope sets.

### Notes

Polynomial zonotopes are in general non-convex. They are always bounded.

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(LazySets.AbstractPolynomialZonotope)
3-element Vector{Any}:
 DensePolynomialZonotope
 SimpleSparsePolynomialZonotope
 SparsePolynomialZonotope
```
"""
abstract type AbstractPolynomialZonotope{N} <: LazySet{N} end

isconvextype(::Type{<:AbstractPolynomialZonotope}) = false
isboundedtype(::Type{<:AbstractPolynomialZonotope}) = true

function Base.:*(M::Union{Real, AbstractMatrix, AbstractVector, LinearAlgebra.UniformScaling}, P::AbstractPolynomialZonotope)
    return linear_map(M, P)
end
