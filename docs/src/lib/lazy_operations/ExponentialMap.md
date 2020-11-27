```@meta
CurrentModule = LazySets
```

# Exponential map

## [Exponential map (ExponentialMap)](@id def_ExponentialMap)

```@docs
ExponentialMap
dim(::ExponentialMap)
ρ(::AbstractVector, ::ExponentialMap)
σ(::AbstractVector, ::ExponentialMap)
∈(::AbstractVector, ::ExponentialMap)
isbounded(::ExponentialMap)
vertices_list(::ExponentialMap)
```
Inherited from [`AbstractAffineMap`](@ref):
* [`an_element`](@ref an_element(::AbstractAffineMap))
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`constraints_list`](@ref constraints_list(::AbstractAffineMap))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractAffineMap))

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

### Sparse matrix exponential

```@docs
SparseMatrixExp
*(::SparseMatrixExp, ::LazySet)
get_row(::SparseMatrixExp, ::Int)
```

## [Exponential projection map (ExponentialProjectionMap)](@id def_ExponentialProjectionMap)

```@docs
ExponentialProjectionMap
dim(::ExponentialProjectionMap)
σ(::AbstractVector, ::ExponentialProjectionMap)
isbounded(::ExponentialProjectionMap)
```
Inherited from [`AbstractAffineMap`](@ref):
* [`an_element`](@ref an_element(::AbstractAffineMap))
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractAffineMap))
* [`vertices_list`](@ref vertices_list(::AbstractAffineMap))
* [`constraints_list`](@ref constraints_list(::AbstractAffineMap))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractAffineMap))

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

### Projection of a sparse matrix exponential

```@docs
ProjectionSparseMatrixExp
*(::ProjectionSparseMatrixExp, ::LazySet)
```
