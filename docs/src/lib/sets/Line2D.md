```@meta
CurrentModule = LazySets
```

# [Line2D](@id def_Line2D)

```@docs
Line2D
dim(::Line2D)
σ(::AbstractVector, ::Line2D)
∈(::AbstractVector, ::Line2D)
an_element(::Line2D{N}) where {N}
rand(::Type{Line2D})
isbounded(::Line2D)
isuniversal(::Line2D, ::Bool=false)
isempty(::Line2D)
constrained_dimensions(::Line2D)
constraints_list(::Line2D)
translate(::Line2D, ::AbstractVector)
project(::AbstractVector, ::Line2D)
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
