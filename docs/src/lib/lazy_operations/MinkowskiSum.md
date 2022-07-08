```@meta
CurrentModule = LazySets
```

# Minkowski sum

## [Binary Minkowski sum (MinkowskiSum)](@id def_MinkowskiSum)

```@docs
MinkowskiSum
⊕(::ConvexSet, ::ConvexSet)
+(::ConvexSet, ::ConvexSet)
swap(::MinkowskiSum)
dim(::MinkowskiSum)
ρ(::AbstractVector, ::MinkowskiSum)
σ(::AbstractVector, ::MinkowskiSum)
isbounded(::MinkowskiSum)
isempty(::MinkowskiSum)
center(::MinkowskiSum)
constraints_list(::MinkowskiSum)
∈(::AbstractVector, ::MinkowskiSum{N, S1, S2}) where {N, S1<:AbstractSingleton, S2<:ConvexSet}
vertices_list(::MinkowskiSum)
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`an_element`](@ref an_element(::ConvexSet{N}) where {N})
* [`singleton_list`](@ref singleton_list(::ConvexSet))

## [``n``-ary Minkowski sum (MinkowskiSumArray)](@id def_MinkowskiSumArray)

```@docs
MinkowskiSumArray
dim(::MinkowskiSumArray)
ρ(::AbstractVector, ::MinkowskiSumArray)
σ(::AbstractVector, ::MinkowskiSumArray)
isbounded(::MinkowskiSumArray)
isempty(::MinkowskiSumArray)
array(::MinkowskiSumArray)
center(::MinkowskiSumArray)
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`an_element`](@ref an_element(::ConvexSet{N}) where {N})
* [`singleton_list`](@ref singleton_list(::ConvexSet))

## [``n``-ary Minkowski sum with cache (CachedMinkowskiSumArray)](@id def_CachedMinkowskiSumArray)

```@docs
CachedMinkowskiSumArray
dim(::CachedMinkowskiSumArray)
σ(::AbstractVector, ::CachedMinkowskiSumArray)
isbounded(::CachedMinkowskiSumArray)
isempty(::CachedMinkowskiSumArray)
array(::CachedMinkowskiSumArray)
forget_sets!(::CachedMinkowskiSumArray)
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`an_element`](@ref an_element(::ConvexSet))
* [`singleton_list`](@ref singleton_list(::ConvexSet))
