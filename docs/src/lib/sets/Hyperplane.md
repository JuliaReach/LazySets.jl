```@meta
CurrentModule = LazySets.HyperplaneModule
```

# [Hyperplane](@id def_Hyperplane)

```@docs
Hyperplane
```

## Operations

```@docs
an_element(::Hyperplane)
constrained_dimensions(::Hyperplane)
constraints_list(::Hyperplane)
dim(::Hyperplane)
isbounded(::Hyperplane)
isempty(::Hyperplane)
isuniversal(::Hyperplane, ::Bool=false)
normalize(::Hyperplane{N}, p::Real=N(2)) where {N}
rand(::Type{Hyperplane})
distance(::AbstractVector, ::Hyperplane{N}) where {N}
∈(::AbstractVector, ::Hyperplane)
project(::AbstractVector, ::Hyperplane)
reflect(::AbstractVector, ::Hyperplane)
ρ(::AbstractVector, ::Hyperplane)
σ(::AbstractVector, ::Hyperplane)
translate(::Hyperplane, ::AbstractVector)
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
