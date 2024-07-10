```@meta
CurrentModule = LazySets.LineModule
```

# [Line](@id def_Line)

```@docs
Line
```

## Operations

```@docs
an_element(::Line)
constraints_list(::Line)
dim(::Line)
direction(::Line)
isbounded(::Line)
isempty(::Line)
isuniversal(::Line; ::Bool=false)
normalize(::Line{N}, ::Real=N(2)) where{N}
normalize!(::Line{N}, ::Real=N(2)) where{N}
rand(::Type{Line})
distance(::AbstractVector, ::Line; ::Real=2.0)
∈(::AbstractVector, ::Line)
linear_map(::AbstractMatrix, ::Line)
ρ(::AbstractVector, ::Line)
σ(::AbstractVector, ::Line)
translate!(::Line, ::AbstractVector)
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
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
