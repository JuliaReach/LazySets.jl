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
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`an_element`](@ref an_element(::ConvexSet{N}) where {N})
* [`singleton_list`](@ref singleton_list(::ConvexSet))

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
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`an_element`](@ref an_element(::ConvexSet{N}) where {N})
* [`singleton_list`](@ref singleton_list(::ConvexSet))
