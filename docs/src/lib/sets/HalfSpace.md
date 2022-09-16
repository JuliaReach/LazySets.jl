```@meta
CurrentModule = LazySets
```

# [Half-space (HalfSpace)](@id def_HalfSpace)

```@docs
HalfSpace
LinearConstraint
dim(::HalfSpace)
ρ(::AbstractVector, ::HalfSpace)
σ(::AbstractVector, ::HalfSpace)
∈(::AbstractVector, ::HalfSpace)
an_element(::HalfSpace)
rand(::Type{HalfSpace})
normalize(::HalfSpace{N}, p=N(2)) where {N}
isbounded(::HalfSpace)
isuniversal(::HalfSpace, ::Bool=false)
isempty(::HalfSpace)
constraints_list(::HalfSpace)
constrained_dimensions(::HalfSpace)
translate(::HalfSpace, ::AbstractVector)
halfspace_left(::AbstractVector, ::AbstractVector)
halfspace_right(::AbstractVector, ::AbstractVector)
complement(::HalfSpace)
iscomplement(::HalfSpace{N}, ::HalfSpace) where {N}
project(::HalfSpace{N}, ::AbstractVector{Int}) where {N}
distance(::AbstractVector, ::HalfSpace{N}) where {N}
```
Inherited from [`LazySet`](@ref):
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))

Inherited from [`AbstractPolyhedron`](@ref):
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolyhedron))
