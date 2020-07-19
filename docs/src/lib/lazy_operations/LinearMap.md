```@meta
CurrentModule = LazySets
```

# [Linear map (LinearMap)](@id def_LinearMap)

```@docs
LinearMap
*(::Union{AbstractMatrix, UniformScaling, AbstractVector, Real}, ::LazySet)
dim(::LinearMap)
ρ(::AbstractVector{N}, ::LinearMap{N}) where {N<:Real}
σ(::AbstractVector{N}, ::LinearMap{N}) where {N<:Real}
∈(::AbstractVector{N}, ::LinearMap{N}) where {N<:Real}
an_element(::LinearMap{N}) where {N<:Real}
vertices_list(::LinearMap{N}) where {N<:Real}
constraints_list(::LinearMap{N}) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::LinearMap{N}) where {N<:Real}
```
Inherited from [`AbstractAffineMap`](@ref):
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

The lazy projection of a set can be conveniently constructed using `Projection`.

```@docs
Projection
```
