```@meta
CurrentModule = LazySets
```

# [Line](@id def_Line)

```@docs
Line
dim(::Line)
ρ(::AbstractVector, ::Line)
σ(::AbstractVector, ::Line)
∈(::AbstractVector, ::Line)
an_element(::Line)
direction(::Line)
rand(::Type{Line})
isbounded(::Line)
isuniversal(::Line; ::Bool=false)
isempty(::Line)
constraints_list(::Line{N, VN}) where {N, VN}
translate(::Line, ::AbstractVector)
translate!(::Line, ::AbstractVector)
normalize(::Line, ::Real=2.0)
normalize!(::Line, ::Real=2.0)
distance(::AbstractVector, ::Line, ::Real=2.0)
distance(::AbstractSingleton, ::Line, ::Real=2.0)
linear_map(::AbstractMatrix, ::Line)
constrained_dimensions(::Line)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
