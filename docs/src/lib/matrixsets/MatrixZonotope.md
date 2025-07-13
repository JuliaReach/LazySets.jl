```@meta
CurrentModule = LazySets.MatrixZonotopeModule
```

# [MatrixZonotope](@id def_MatrixZonotope)

```@docs
MatrixZonotope
```

## Operations

```@docs
center(::MatrixZonotope)
generators(::MatrixZonotope)
transpose(::MatrixZonotope)
size(::MatrixZonotope)
ngens(::MatrixZonotope)
rand(::Type{MatrixZonotope})
norm(::MatrixZonotope, ::Real)
_rowwise_zonotope_norm
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`eltype`](@ref eltype(::Type{<:LazySet}))