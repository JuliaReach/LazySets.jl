```@meta
CurrentModule = LazySets
```

# [Line2D](@id def_Line2D)

```@docs
Line2D
dim(::Line2D)
σ(::AbstractVector, ::Line2D)
∈(::AbstractVector, ::Line2D)
an_element(::Line2D)
rand(::Type{Line2D})
isbounded(::Line2D)
isuniversal(::Line2D, ::Bool=false)
isempty(::Line2D)
constrained_dimensions(::Line2D)
constraints_list(::Line2D)
translate(::Line2D, ::AbstractVector)
project(::AbstractVector, ::Line2D)
```
Inherited from [`LazySet`](@ref):
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`low`](@ref low(::LazySet, ::Int))
* [`high`](@ref high(::LazySet, ::Int))
* [`reflect`](@ref reflect(::LazySet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolyhedron))
