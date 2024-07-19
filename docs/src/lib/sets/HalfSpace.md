```@meta
CurrentModule = LazySets.HalfSpaceModule
```

# [Half-space (HalfSpace)](@id def_HalfSpace)

```@docs
HalfSpace
LinearConstraint
```

## Operations

```@docs
an_element(::HalfSpace)
complement(::HalfSpace)
constrained_dimensions(::HalfSpace)
constraints_list(::HalfSpace)
dim(::HalfSpace)
isbounded(::HalfSpace)
isempty(::HalfSpace)
isuniversal(::HalfSpace, ::Bool=false)
normalize(::HalfSpace{N}, p::Real=N(2)) where {N}
rand(::Type{HalfSpace})
distance(::AbstractVector, ::HalfSpace{N}) where {N}
∈(::AbstractVector, ::HalfSpace)
halfspace_left(::AbstractVector, ::AbstractVector)
halfspace_right(::AbstractVector, ::AbstractVector)
iscomplement(::HalfSpace{N}, ::HalfSpace) where {N}
project(::HalfSpace{N}, ::AbstractVector{Int}) where {N}
ρ(::AbstractVector, ::HalfSpace)
σ(::AbstractVector, ::HalfSpace)
translate(::HalfSpace, ::AbstractVector)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`reflect`](@ref reflect(::LazySet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolyhedron))
