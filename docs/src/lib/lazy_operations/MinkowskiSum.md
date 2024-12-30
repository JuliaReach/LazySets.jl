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
center(::MinkowskiSum)
constraints_list(::MinkowskiSum)
∈(::AbstractVector, ::MinkowskiSum{N,<:AbstractSingleton}) where {N}
vertices_list(::MinkowskiSum)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet)
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`reflect`](@ref reflect(::LazySet))

## [``n``-ary Minkowski sum (MinkowskiSumArray)](@id def_MinkowskiSumArray)

```@docs
MinkowskiSumArray
⊕(::LazySet, ::LazySet...)
+(::LazySet, ::LazySet...)
dim(::MinkowskiSumArray)
ρ(::AbstractVector, ::MinkowskiSumArray)
σ(::AbstractVector, ::MinkowskiSumArray)
isbounded(::MinkowskiSumArray)
isempty(::MinkowskiSumArray)
array(::MinkowskiSumArray)
center(::MinkowskiSumArray)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet)
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`reflect`](@ref reflect(::LazySet))

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
* [`reflect`](@ref reflect(::LazySet))
