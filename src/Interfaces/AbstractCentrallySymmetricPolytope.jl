export AbstractCentrallySymmetricPolytope

"""
    AbstractCentrallySymmetricPolytope{N} <: AbstractPolytope{N}

Abstract type for centrally symmetric, polytopic sets.
It combines the `AbstractCentrallySymmetric` and `AbstractPolytope` interfaces.
Such a type combination is necessary as long as Julia does not support
[multiple inheritance](https://github.com/JuliaLang/julia/issues/5).

### Notes

Every concrete `AbstractCentrallySymmetricPolytope` must define the following
`AbstractCentrallySymmetric` method, in addition to the `AbstractPolytope`
methods:

- `center(::AbstractCentrallySymmetricPolytope)` -- return the center point

The subtypes of `AbstractCentrallySymmetricPolytope` (including abstract
interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractCentrallySymmetricPolytope)
2-element Vector{Any}:
 AbstractZonotope
 Ball1
```
"""
abstract type AbstractCentrallySymmetricPolytope{N} <: AbstractPolytope{N} end
