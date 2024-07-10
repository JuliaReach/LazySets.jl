export AbstractPolynomialZonotope

"""
    AbstractPolynomialZonotope{N} <: LazySet{N}

Abstract type for polynomial zonotope sets.

### Notes

Polynomial zonotopes are in general non-convex. They are always bounded.

Every concrete `AbstractPolynomialZonotope` must define the following functions:

- `center(::AbstractPolynomialZonotope)` -- return the center

- `order(::AbstractPolynomialZonotope)` -- return the order

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractPolynomialZonotope)
3-element Vector{Any}:
 DensePolynomialZonotope
 SimpleSparsePolynomialZonotope
 SparsePolynomialZonotope
```
"""
abstract type AbstractPolynomialZonotope{N} <: LazySet{N} end

"""
    order(PZ::AbstractPolynomialZonotope)

Return the order of a polynomial zonotope.

### Input

- `PZ` -- polynomial zonotope

### Output

A rational number representing the order of `PZ`.

### Notes

The order of a polynomial zonotope is defined as the quotient of its number of
generators and its dimension.
"""
function order(::AbstractPolynomialZonotope) end

isconvextype(::Type{<:AbstractPolynomialZonotope}) = false
isboundedtype(::Type{<:AbstractPolynomialZonotope}) = true

"""
    dim(PZ::AbstractPolynomialZonotope)

Return the ambient dimension of a polynomial zonotope.

### Input

- `PZ` -- polynomial zonotope

### Output

An integer representing the ambient dimension of the polynomial zonotope.
"""
dim(PZ::AbstractPolynomialZonotope) = length(center(PZ))
