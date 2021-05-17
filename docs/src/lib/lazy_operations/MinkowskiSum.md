```@meta
CurrentModule = LazySets
```

# Minkowski sum

## [Binary Minkowski sum (MinkowskiSum)](@id def_MinkowskiSum)

```@docs
MinkowskiSum
⊕(::LazySet, ::LazySet)
+(::LazySet, ::LazySet)
swap(::MinkowskiSum)
dim(::MinkowskiSum)
ρ(::AbstractVector, ::MinkowskiSum)
σ(::AbstractVector, ::MinkowskiSum)
isbounded(::MinkowskiSum)
isempty(::MinkowskiSum)
constraints_list(::MinkowskiSum)
∈(::AbstractVector, ::MinkowskiSum{N, S1, S2}) where {N, S1<:AbstractSingleton, S2<:LazySet}
vertices_list(::MinkowskiSum)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N})
* [`singleton_list`](@ref singleton_list(::LazySet))

## [``n``-ary Minkowski sum (MinkowskiSumArray)](@id def_MinkowskiSumArray)

```@docs
MinkowskiSumArray
dim(::MinkowskiSumArray)
ρ(::AbstractVector, ::MinkowskiSumArray)
σ(::AbstractVector, ::MinkowskiSumArray)
isbounded(::MinkowskiSumArray)
isempty(::MinkowskiSumArray)
array(::MinkowskiSumArray)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N})
* [`singleton_list`](@ref singleton_list(::LazySet))

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
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
