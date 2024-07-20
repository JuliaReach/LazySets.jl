```@meta
CurrentModule = LazySets
```

# [Line2D](@id def_Line2D)

```@docs
Line2D
```

## Operations

```@docs
an_element(::Line2D)
constrained_dimensions(::Line2D)
constraints_list(::Line2D)
dim(::Line2D)
isbounded(::Line2D)
isempty(::Line2D)
isuniversal(::Line2D, ::Bool=false)
rand(::Type{Line2D})
∈(::AbstractVector, ::Line2D)
project(::AbstractVector, ::Line2D)
σ(::AbstractVector, ::Line2D)
translate(::Line2D, ::AbstractVector)
intersection(::Line2D, ::Line2D)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`reflect`](@ref reflect(::LazySet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolyhedron))
