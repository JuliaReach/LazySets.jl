```@contents
Pages = ["template_directions.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets.Approximations
```

# Template directions

See also [`overapproximate(X::LazySet, dir::AbstractDirections)`](@ref).

```@docs
AbstractDirections
isbounding
isnormalized
project(::LazySet, ::AbstractVector{Int}, ::Type{<:AbstractDirections})
BoxDirections
DiagDirections
OctDirections
BoxDiagDirections
PolarDirections
SphericalDirections
CustomDirections
```
