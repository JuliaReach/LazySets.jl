```@meta
CurrentModule = LazySets.DensePolynomialZonotopeModule
```

# [Polynomial zonotope (DensePolynomialZonotope)](@id def_DensePolynomialZonotope)

```@docs
DensePolynomialZonotope
```

## Operations

```@docs
center(::DensePolynomialZonotope)
ngens_dep(::DensePolynomialZonotope)
ngens_indep(::DensePolynomialZonotope)
polynomial_order(::DensePolynomialZonotope)
linear_map(::AbstractMatrix, ::DensePolynomialZonotope)
scale!(::Real, ::DensePolynomialZonotope)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`AbstractPolynomialZonotope`](@ref):
* [`dim`](@ref dim(::AbstractPolynomialZonotope))
* [`order`](@ref dim(::AbstractPolynomialZonotope))
