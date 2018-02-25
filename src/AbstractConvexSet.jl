export AbstractConvexSet,
       isconvex

"""
    AbstractConvexSet{N} <: LazySet{N}

Abstract type for convex sets, i.e., sets characterized by a (possibly infinite)
intersection of halfspaces, or equivalently, sets ``S`` such that for any two
elements ``x, y ∈ S`` and ``0 ≤ λ ≤ 1`` it holds that ``λ x + (1-λ) y ∈ S``.

### Notes

```jldoctest
julia> subtypes(AbstractConvexSet)
6-element Array{Union{DataType, UnionAll},1}:
 LazySets.AbstractPointSymmetric
 LazySets.AbstractPolytope
 LazySets.EmptySet
 LazySets.HalfSpace
 LazySets.Hyperplane
 LazySets.Line
```
"""
abstract type AbstractConvexSet{N} <: LazySet{N} end


"""
    isconvex(S::AbstractConvexSet)::Bool

Sufficient check if a set is convex.

### Input

- `S` -- convex set

### Output

`true`.
"""
function isconvex(S::AbstractConvexSet)::Bool
    return true
end
