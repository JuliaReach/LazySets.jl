```@meta
CurrentModule = LazySets
```

# [Line](@id def_Line)

```@docs
Line
dim(::Line)
σ(::AbstractVector{N}, ::Line{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Line{N}) where {N<:Real}
an_element(::Line{N}) where {N<:Real}
rand(::Type{Line})
isbounded(::Line)
isuniversal(::Line{N}, ::Bool=false) where {N<:Real}
isempty(::Line)
constrained_dimensions(::Line{N}) where {N<:Real}
constraints_list(::Line{N}) where {N<:Real}
translate(::Line{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
