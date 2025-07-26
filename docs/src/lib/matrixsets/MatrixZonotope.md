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
```

# [MatrixZonotope](@id def_MatrixZonotope)

```@docs
MatrixZonotope
```

## Operations

```@docs
center(::MatrixZonotope)
generators(::MatrixZonotope)
indexvector(::MatrixZonotope)
transpose(::MatrixZonotope)
ngens(::MatrixZonotope)
rand(::Type{MatrixZonotope})
*(::Real, ::MatrixZonotope)
norm(::MatrixZonotope, ::Real)
linear_map(::AbstractMatrix, ::MatrixZonotope)
linear_map(::MatrixZonotope, ::AbstractMatrix)
_rowwise_zonotope_norm
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
factors
nfactors
remove_redundant_factors
```

```@meta
CurrentModule = LazySets.MatrixZonotopeModule
```

## Operations 
Inherited from [`AbstractMatrixZonotope`](@ref):
* [`size`](@ref size(::AbstractMatrixZonotope))

# [MatrixZonotopeExp](@id def_MatrixZonotopeExp)
```@docs
MatrixZonotopeExp
```
## Operations 
Inherited from [`AbstractMatrixZonotope`](@ref):
* [`size`](@ref size(::AbstractMatrixZonotope))