```@meta
CurrentModule = LazySets
```

# [Inverse linear map (InverseLinearMap)](@id def_InverseLinearMap)

```@docs
InverseLinearMap
dim(::InverseLinearMap)
σ(::AbstractVector{N}, ::InverseLinearMap{N,S,NM,MAT}) where {N<:Real, S<:LazySet{N}, NM, MAT}
ρ(::AbstractVector{N}, ::InverseLinearMap{N,S,NM,MAT}) where {N<:Real, S<:LazySet{N}, NM, MAT}
∈(::AbstractVector, ::InverseLinearMap)
an_element(::InverseLinearMap)
vertices_list(::InverseLinearMap{N})  where {N}
constraints_list(::InverseLinearMap)
linear_map(::AbstractMatrix{N}, ::InverseLinearMap{N}) where {N<:Real}
```
Inherited from [`AbstractAffineMap`](@ref):
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))
