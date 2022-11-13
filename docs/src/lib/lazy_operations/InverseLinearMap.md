```@meta
CurrentModule = LazySets
```

# [Inverse linear map (InverseLinearMap)](@id def_InverseLinearMap)

```@docs
InverseLinearMap
dim(::InverseLinearMap)
σ(::AbstractVector, ::InverseLinearMap)
ρ(::AbstractVector, ::InverseLinearMap)
∈(::AbstractVector, ::InverseLinearMap)
an_element(::InverseLinearMap)
vertices_list(::InverseLinearMap)
constraints_list(::InverseLinearMap)
linear_map(::AbstractMatrix, ::InverseLinearMap)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractAffineMap`](@ref):
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))
