```@meta
CurrentModule = LazySets
```

# Cartesian product

## [Binary Cartesian product (CartesianProduct)](@id def_CartesianProduct)

```@docs
CartesianProduct
×(::LazySet, ::LazySet)
*(::LazySet, ::LazySet)
swap(::CartesianProduct)
dim(::CartesianProduct)
ρ(::AbstractVector, ::CartesianProduct)
σ(::AbstractVector, ::CartesianProduct)
isbounded(::CartesianProduct)
∈(::AbstractVector, ::CartesianProduct)
isempty(::CartesianProduct)
constraints_list(::CartesianProduct)
vertices_list(::CartesianProduct)
linear_map(M::AbstractMatrix, cp::CartesianProduct)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet)
* [`singleton_list`](@ref singleton_list(::LazySet))

## [``n``-ary Cartesian product (CartesianProductArray)](@id def_CartesianProductArray)

```@docs
CartesianProductArray
dim(::CartesianProductArray)
ρ(::AbstractVector, ::CartesianProductArray)
σ(::AbstractVector, ::CartesianProductArray)
isbounded(::CartesianProductArray)
∈(::AbstractVector, ::CartesianProductArray)
isempty(::CartesianProductArray)
constraints_list(::CartesianProductArray)
vertices_list(::CartesianProductArray)
linear_map(M::AbstractMatrix, cpa::CartesianProductArray)
array(::CartesianProductArray{)
block_structure(cpa::CartesianProductArray)
block_to_dimension_indices(cpa::CartesianProductArray, vars::Vector{Int})
substitute_blocks(low_dim_cpa::CartesianProductArray, orig_cpa::CartesianProductArray,
blocks::Vector{Tuple{Int,Int}})
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N})
* [`singleton_list`](@ref singleton_list(::LazySet))
