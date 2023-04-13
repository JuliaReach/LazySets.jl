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
constraints_list(::Line)
translate(::Line, ::AbstractVector)
translate!(::Line, ::AbstractVector)
normalize(::Line{N}, ::Real=N(2)) where{N}
normalize!(::Line{N}, ::Real=N(2)) where{N}
distance(::AbstractVector, ::Line; ::Real=2.0)
linear_map(::AbstractMatrix, ::Line)
```
Inherited from [`LazySet`](@ref):
* [`low`](@ref low(::LazySet))
* [`high`](@ref high(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`low`](@ref low(::LazySet, ::Int))
* [`high`](@ref high(::LazySet, ::Int))
* [`reflect`](@ref reflect(::LazySet))
