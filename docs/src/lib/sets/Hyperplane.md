```@meta
CurrentModule = LazySets
```

# [Hyperplane](@id def_Hyperplane)

```@docs
Hyperplane
dim(::Hyperplane)
ρ(::AbstractVector, ::Hyperplane)
σ(::AbstractVector, ::Hyperplane)
∈(::AbstractVector, ::Hyperplane)
an_element(::Hyperplane)
rand(::Type{Hyperplane})
isbounded(::Hyperplane)
isuniversal(::Hyperplane, ::Bool=false)
isempty(::Hyperplane)
constrained_dimensions(::Hyperplane)
constraints_list(::Hyperplane)
translate(::Hyperplane, ::AbstractVector)
normalize(::Hyperplane{N}, p::Real=N(2)) where {N}
distance(::AbstractVector, ::Hyperplane{N}) where {N}
reflect(::AbstractVector, ::Hyperplane)
project(::AbstractVector, ::Hyperplane)
```
Inherited from [`LazySet`](@ref):
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`reflect`](@ref reflect(::LazySet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolyhedron))
