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
ngens(::MatrixZonotope)
norm(::MatrixZonotope, ::Real)
_rowwise_zonotope_norm
rand(::Type{MatrixZonotope})
size(::MatrixZonotope)
transpose(::MatrixZonotope)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
