```@meta
CurrentModule = LazySets
```

# Convex hull

## [Binary convex hull (ConvexHull)](@id def_ConvexHull)

```@docs
ConvexHull
CH
swap(::ConvexHull)
dim(::ConvexHull)
ρ(::AbstractVector{N}, ::ConvexHull{N}) where {N<:Real}
σ(::AbstractVector{N}, ::ConvexHull{N}) where {N<:Real}
isbounded(::ConvexHull)
isempty(::ConvexHull)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

## [``n``-ary convex hull (ConvexHullArray)](@id def_ConvexHullArray)

```@docs
ConvexHullArray
CHArray
dim(::ConvexHullArray)
ρ(::AbstractVector{N}, ::ConvexHullArray{N}) where {N<:Real}
σ(::AbstractVector{N}, ::ConvexHullArray{N}) where {N<:Real}
isbounded(::ConvexHullArray)
array(::ConvexHullArray{N, S}) where {N<:Real, S<:LazySet{N}}
isempty(::ConvexHullArray)
vertices_list(::ConvexHullArray{N, Singleton{N, VT}}) where {N, VT}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})
