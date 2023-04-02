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
ρ(::AbstractVector, ::ConvexHull)
σ(::AbstractVector, ::ConvexHull)
isbounded(::ConvexHull)
isempty(::ConvexHull)
vertices_list(::ConvexHull)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet)
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`reflect`](@ref reflect(::LazySet))

## [``n``-ary convex hull (ConvexHullArray)](@id def_ConvexHullArray)

```@docs
ConvexHullArray
CHArray
dim(::ConvexHullArray)
ρ(::AbstractVector, ::ConvexHullArray)
σ(::AbstractVector, ::ConvexHullArray)
isbounded(::ConvexHullArray)
array(::ConvexHullArray)
isempty(::ConvexHullArray)
vertices_list(::ConvexHullArray)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet)
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
