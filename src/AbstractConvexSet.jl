export AbstractConvexSet, LazySet

"""
    AbstractConvexSet{N}

Abstract type for convex sets, i.e., sets characterized by a (possibly infinite)
intersection of halfspaces, or equivalently, sets ``S`` such that for any two
elements ``x, y ∈ S`` and ``0 ≤ λ ≤ 1`` it holds that ``λ x + (1-λ) y ∈ S``.

### Notes

`AbstractConvexSet` types should be parameterized with a type `N`, typically
`N<:Real`, for using different numeric types.

Every concrete `AbstractConvexSet` must define the following functions:
- `σ(d::AbstractVector{N}, S::AbstractConvexSet)::AbstractVector{N}` -- the
    support vector of `S` in a given direction `d`
- `dim(S::AbstractConvexSet)::Int` -- the ambient dimension of `S`

```jldoctest
julia> subtypes(AbstractConvexSet)
15-element Array{Union{DataType, UnionAll},1}:
 LazySets.AbstractPointSymmetric  
 LazySets.AbstractPolytope        
 LazySets.CartesianProduct        
 LazySets.CartesianProductArray   
 LazySets.ConvexHull              
 LazySets.ConvexHullArray         
 LazySets.EmptySet                
 LazySets.ExponentialMap          
 LazySets.ExponentialProjectionMap
 LazySets.HalfSpace               
 LazySets.Hyperplane              
 LazySets.Intersection            
 LazySets.LinearMap               
 LazySets.MinkowskiSum            
 LazySets.MinkowskiSumArray
```
"""
abstract type AbstractConvexSet{N} end

"""
    LazySet

Alias for `AbstractConvexSet`.
"""
const LazySet = AbstractConvexSet
