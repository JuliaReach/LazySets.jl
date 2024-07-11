```@meta
CurrentModule = LazySets
```

# [Reset map (ResetMap)](@id def_ResetMap)

```@docs
ResetMap
dim(::ResetMap)
ρ(::AbstractVector, ::ResetMap)
σ(::AbstractVector, ::ResetMap)
an_element(::ResetMap)
matrix(::ResetMap{N}) where {N}
vector(::ResetMap)
set(::ResetMap)
constraints_list(::ResetMap)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`reflect`](@ref reflect(::LazySet))

Inherited from [`AbstractAffineMap`](@ref):
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractAffineMap))
* [`vertices_list`](@ref vertices_list(::AbstractAffineMap))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractAffineMap))
