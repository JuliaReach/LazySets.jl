```@meta
CurrentModule = LazySets
```

# Exponential map

## [Exponential map (ExponentialMap)](@id def_ExponentialMap)

```@docs
ExponentialMap
dim(::ExponentialMap)
ρ(::AbstractVector{N}, ::ExponentialMap{N}) where {N<:Real}
σ(::AbstractVector{N}, ::ExponentialMap{N}) where {N<:Real}
∈(::AbstractVector{N}, ::ExponentialMap{N}) where {N<:Real}
isbounded(::ExponentialMap)
vertices_list(::ExponentialMap{N}) where {N<:Real}
```
Inherited from [`AbstractAffineMap`](@ref):
* [`an_element`](@ref an_element(::AbstractAffineMap))
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`constraints_list`](@ref constraints_list(::AbstractAffineMap{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractAffineMap{N}) where {N<:Real})

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

### Sparse matrix exponential

```@docs
SparseMatrixExp
*(::SparseMatrixExp{N}, ::LazySet{N}) where {N<:Real}
get_row(::SparseMatrixExp, ::Int)
```

## [Exponential projection map (ExponentialProjectionMap)](@id def_ExponentialProjectionMap)

```@docs
ExponentialProjectionMap
dim(::ExponentialProjectionMap)
σ(::AbstractVector{N}, ::ExponentialProjectionMap{N}) where {N<:Real}
isbounded(::ExponentialProjectionMap)
```
Inherited from [`AbstractAffineMap`](@ref):
* [`an_element`](@ref an_element(::AbstractAffineMap))
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractAffineMap{N}) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractAffineMap{N}) where {N<:Real})
* [`constraints_list`](@ref constraints_list(::AbstractAffineMap{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractAffineMap{N}) where {N<:Real})

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

### Projection of a sparse matrix exponential

```@docs
ProjectionSparseMatrixExp
*(::ProjectionSparseMatrixExp, ::LazySet)
```
