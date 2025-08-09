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
order(::MatrixZonotope)
transpose(::MatrixZonotope)
ngens(::MatrixZonotope)
rand(::Type{MatrixZonotope})
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

```@meta
CurrentModule = LazySets.MatrixZonotopeModule
```

Undocumented implementations:
* [`size`](@ref size(::AbstractMatrixZonotope))
* [`remove_redundant_generators`](@ref remove_redundant_generators(::AbstractZonotope))
* [`minkowski_sum`](@ref minkowski_sum(::LazySet, LazySet))
* [`reduce_order`](@ref reduce_order(::AbstractZonotope, ::Real, ::AbstractReductionMethod=GIR05()))

# [MatrixZonotopeProduct](@id def_MatrixZonotopeProduct)
```@docs
MatrixZonotopeProduct
factors
nfactors
```

## Operations
Undocumented implementations:
* [`size`](@ref size(::AbstractMatrixZonotope))

# [MatrixZonotopeExp](@id def_MatrixZonotopeExp)
```@docs
MatrixZonotopeExp
```
## Operations
Undocumented implementations:
* [`size`](@ref size(::AbstractMatrixZonotope))
