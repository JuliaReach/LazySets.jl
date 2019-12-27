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
ρ(::AbstractVector{N}, ::CartesianProduct{N}) where {N<:Real}
σ(::AbstractVector{N}, ::CartesianProduct{N}) where {N<:Real}
isbounded(::CartesianProduct)
∈(::AbstractVector{N}, ::CartesianProduct{N}) where {N<:Real}
isempty(::CartesianProduct)
constraints_list(::CartesianProduct{N}) where {N<:Real}
vertices_list(::CartesianProduct{N}) where {N<:Real}
linear_map(M::AbstractMatrix{N}, cp::CartesianProduct{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

## [``n``-ary Cartesian product (CartesianProductArray)](@id def_CartesianProductArray)

```@docs
CartesianProductArray
dim(::CartesianProductArray)
ρ(::AbstractVector{N}, ::CartesianProductArray{N}) where {N<:Real}
σ(::AbstractVector{N}, ::CartesianProductArray{N}) where {N<:Real}
isbounded(::CartesianProductArray)
∈(::AbstractVector{N}, ::CartesianProductArray{N}) where {N<:Real}
isempty(::CartesianProductArray)
constraints_list(::CartesianProductArray{N}) where {N<:Real}
vertices_list(::CartesianProductArray{N}) where {N<:Real}
linear_map(M::AbstractMatrix{N}, cpa::CartesianProductArray{N}) where {N<:Real}
array(::CartesianProductArray{N, S}) where {N<:Real, S<:LazySet{N}}
block_structure(cpa::CartesianProductArray{N}) where {N}
block_to_dimension_indices(cpa::CartesianProductArray{N}, vars::Vector{Int}) where {N}
substitute_blocks(low_dim_cpa::CartesianProductArray{N}, orig_cpa::CartesianProductArray{N},
blocks::Vector{Tuple{Int,Int}}) where {N}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})
