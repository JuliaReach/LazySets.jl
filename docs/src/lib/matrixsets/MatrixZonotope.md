```@meta
CurrentModule = LazySets.MatrixZonotopeModule
```

# [AbstractMatrixZonotope](@id def_MatrixZonotopeProduct)
```@docs
AbstractMatrixZonotope
```

This interface requires to implement the following functions:
```@docs
size(::AbstractMatrixZonotope)
size(::AbstractMatrixZonotope, ::Int)
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
ngens(::MatrixZonotope)
rand(::Type{MatrixZonotope})
*(::Real, ::MatrixZonotope)
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
* [`eltype`](@ref eltype(::Type{<:LazySet}))

```@meta
CurrentModule = LazySets.MatrixZonotopeModule
```

Inherited from [`AbstractMatrixZonotope`](@ref):
* [`size`](@ref size(::AbstractMatrixZonotope))

```@meta
CurrentModule = LazySets.MatrixZonotopeModule
```

# [MatrixZonotopeProduct](@id def_MatrixZonotopeProduct)
```@docs
MatrixZonotopeProduct
*(::MatrixZonotope, ::MatrixZonotope)
```

```@meta
CurrentModule = LazySets.MatrixZonotopeModule
```

## Operations 
Inherited from [`AbstractMatrixZonotope`](@ref):
* [`size`](@ref size(::AbstractMatrixZonotope))
