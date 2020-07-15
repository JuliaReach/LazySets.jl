```@meta
CurrentModule = LazySets
```

# [Line2D](@id def_Line2D)

```@docs
Line2D
dim(::Line2D)
σ(::AbstractVector{N}, ::Line2D{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Line2D{N}) where {N<:Real}
an_element(::Line2D{N}) where {N<:Real}
rand(::Type{Line2D})
isbounded(::Line2D)
isuniversal(::Line2D{N}, ::Bool=false) where {N<:Real}
isempty(::Line2D)
constrained_dimensions(::Line2D{N}) where {N<:Real}
constraints_list(::Line2D{N}) where {N<:Real}
translate(::Line2D{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
