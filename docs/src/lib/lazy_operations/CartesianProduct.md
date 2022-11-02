```@meta
CurrentModule = LazySets
```

# Cartesian product

## [Binary Cartesian product (CartesianProduct)](@id def_CartesianProduct)

```@docs
CartesianProduct
×(::ConvexSet, ::ConvexSet)
*(::ConvexSet, ::ConvexSet)
swap(::CartesianProduct)
dim(::CartesianProduct)
ρ(::AbstractVector, ::CartesianProduct)
σ(::AbstractVector, ::CartesianProduct)
isbounded(::CartesianProduct)
∈(::AbstractVector, ::CartesianProduct)
isempty(::CartesianProduct)
center(::CartesianProduct)
constraints_list(::CartesianProduct)
vertices_list(::CartesianProduct)
linear_map(::AbstractMatrix, ::CartesianProduct)
volume(::CartesianProduct)
project(::CartesianProduct{N, IT, HT}, ::AbstractVector{Int}) where {N, IT<:Interval, HT<:AbstractHyperrectangle{N}}
project(::CartesianProduct{N, IT, ZT}, ::AbstractVector{Int}) where {N, IT<:Interval, ZT<:AbstractZonotope{N}}
project(::CartesianProduct{N, IT, Union{VP1, VP2}}, ::AbstractVector{Int}) where {N, IT<:Interval, VP1<:VPolygon{N}, VP2<:VPolytope{N}}
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`an_element`](@ref an_element(::ConvexSet)
* [`singleton_list`](@ref singleton_list(::ConvexSet))

## [``n``-ary Cartesian product (CartesianProductArray)](@id def_CartesianProductArray)

```@docs
CartesianProductArray
dim(::CartesianProductArray)
ρ(::AbstractVector, ::CartesianProductArray)
σ(::AbstractVector, ::CartesianProductArray)
isbounded(::CartesianProductArray)
∈(::AbstractVector, ::CartesianProductArray)
isempty(::CartesianProductArray)
center(::CartesianProductArray)
constraints_list(cpa::CartesianProductArray{N}) where {N}
vertices_list(::CartesianProductArray{N}) where {N}
linear_map(M::AbstractMatrix, cpa::CartesianProductArray)
array(::CartesianProductArray)
block_structure(cpa::CartesianProductArray)
block_to_dimension_indices(cpa::CartesianProductArray, vars::Vector{Int})
substitute_blocks(low_dim_cpa::CartesianProductArray{N}, orig_cpa::CartesianProductArray{N}, blocks::Vector{Tuple{Int,Int}}) where{N}
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`an_element`](@ref an_element(::ConvexSet{N}) where {N})
* [`singleton_list`](@ref singleton_list(::ConvexSet))
