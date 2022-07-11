export LazySet

"""
    LazySet{N}

Abstract type for the set types in LazySets.jl

### Notes

`LazySet` types should be parameterized with a type `N`, typically `N<:Real`,
for using different numeric types.

Every concrete `LazySet` must define the following functions:
- `σ(d::AbstractVector, S::LazySet)` -- the support vector of `S` in a given
    direction `d`
- `dim(S::LazySet)` -- the ambient dimension of `S`

The function
- `ρ(d::AbstractVector, S::LazySet)` -- the support function of `S` in a given
    direction `d`
is optional because there is a fallback implementation relying on `σ`.
However, for unbounded sets (which includes most lazy set types) this fallback
cannot be used and an explicit method must be implemented.

The subtypes of `LazySet` (including abstract interfaces):

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(LazySet, false)
1-element Vector{Any}:
 ConvexSet
```

"""
abstract type LazySet{N} end
