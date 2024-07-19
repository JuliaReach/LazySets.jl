export AbstractPolynomialZonotope,
       ngens, ngens_dep, ngens_indep,
       order,
       polynomial_order

"""
    AbstractPolynomialZonotope{N} <: LazySet{N}

Abstract type for polynomial zonotope sets.

### Notes

Polynomial zonotopes are in general non-convex. They are always bounded.

Every concrete `AbstractPolynomialZonotope` must define the following functions:

- `center(::AbstractPolynomialZonotope)` -- return the center

- `polynomial_order(::AbstractPolynomialZonotope)` -- return the polynomial order

- `ngens_dep` -- return the number of dependent generators

- `ngens_indep` -- return the number of independent generators

The subtypes of `AbstractPolynomialZonotope` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractPolynomialZonotope)
2-element Vector{Any}:
 AbstractSparsePolynomialZonotope
 DensePolynomialZonotope
```
"""
abstract type AbstractPolynomialZonotope{N} <: LazySet{N} end

"""
    polynomial_order(P::AbstractPolynomialZonotope)

Determine the polynomial order of a polynomial zonotope.

### Input

- `P` -- polynomial zonotope

### Output

A nonnegative integer representing the polynomial order.

### Notes

The polynomial order is the maximum sum of all monomials' parameter exponents.
"""
function polynomial_order(::AbstractPolynomialZonotope) end

"""
    ngens_dep(P::AbstractPolynomialZonotope)

Determine the number of dependent generators of a polynomial zonotope.

### Input

- `P` -- polynomial zonotope

### Output

A nonnegative integer representing the number of dependent generators.
"""
function ngens_dep(::AbstractPolynomialZonotope) end

"""
    ngens_indep(P::AbstractPolynomialZonotope)

Determine the number of independent generators of a polynomial zonotope.

### Input

- `P` -- polynomial zonotope

### Output

A nonnegative integer representing the number of independent generators.
"""
function ngens_indep(::AbstractPolynomialZonotope) end

"""
    ngens(P::AbstractPolynomialZonotope)

Return the number of generators of a polynomial zonotope.

### Input

- `P` -- polynomial zonotope

### Output

The total number of generators of `P`.
"""
function ngens(P::AbstractPolynomialZonotope)
    return ngens_dep(P) + ngens_indep(P)
end

"""
    order(P::AbstractPolynomialZonotope)

Return the order of a polynomial zonotope.

### Input

- `P` -- polynomial zonotope

### Output

The order, defined as the quotient between the total number of generators and
the ambient dimension, as a `Rational` number.
"""
function order(P::AbstractPolynomialZonotope)
    return ngens(P) // dim(P)
end

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
