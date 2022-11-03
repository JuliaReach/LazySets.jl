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
matrix(::ResetMap)
vector(::ResetMap)
set(::ResetMap)
constraints_list(::ResetMap)
```
Inherited from [`AbstractAffineMap`](@ref):
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractAffineMap))
* [`vertices_list`](@ref vertices_list(::AbstractAffineMap))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractAffineMap))

Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`singleton_list`](@ref singleton_list(::ConvexSet))
