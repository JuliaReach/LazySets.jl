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
ρ(::AbstractVector{N}, ::MinkowskiSum{N}) where {N<:Real}
σ(::AbstractVector{N}, ::MinkowskiSum{N}) where {N<:Real}
isbounded(::MinkowskiSum)
isempty(::MinkowskiSum)
constraints_list(::MinkowskiSum)
∈(::AbstractVector{N}, ::MinkowskiSum{N, S1, S2}) where {N, S1<:AbstractSingleton, S2<:LazySet}
vertices_list(MS::MinkowskiSum{N, Z1, Z2}) where{N<:Real, Z1<:AbstractZonotope{N}, Z2<:AbstractZonotope{N}}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})
* [`singleton_list`](@ref singleton_list(::LazySet))

## [``n``-ary Minkowski sum (MinkowskiSumArray)](@id def_MinkowskiSumArray)

```@docs
MinkowskiSumArray
dim(::MinkowskiSumArray)
ρ(::AbstractVector{N}, ::MinkowskiSumArray{N}) where {N<:Real}
σ(::AbstractVector{N}, ::MinkowskiSumArray{N}) where {N<:Real}
isbounded(::MinkowskiSumArray)
isempty(::MinkowskiSumArray)
array(::MinkowskiSumArray{N, S}) where {N<:Real, S<:LazySet{N}}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})
* [`singleton_list`](@ref singleton_list(::LazySet))

## [``n``-ary Minkowski sum with cache (CachedMinkowskiSumArray)](@id def_CachedMinkowskiSumArray)

```@docs
CachedMinkowskiSumArray
dim(::CachedMinkowskiSumArray)
σ(::AbstractVector{N}, ::CachedMinkowskiSumArray{N}) where {N<:Real}
isbounded(::CachedMinkowskiSumArray)
isempty(::CachedMinkowskiSumArray)
array(::CachedMinkowskiSumArray{N, S}) where {N<:Real, S<:LazySet{N}}
forget_sets!(::CachedMinkowskiSumArray)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})
* [`singleton_list`](@ref singleton_list(::LazySet))
