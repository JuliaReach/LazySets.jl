```@contents
Pages = ["AbstractPolynomialZonotope.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Polynomial zonotope sets (AbstractPolynomialZonotope)](@id def_AbstractPolynomialZonotope)

```@docs
AbstractPolynomialZonotope
```

This interface requires to implement the following functions:

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
center(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
order(::AbstractPolynomialZonotope)
```

This interface defines the following functions:

```@docs
dim(::AbstractPolynomialZonotope)
```

## Implementations

* [Dense polynomial zonotope (DensePolynomialZonotope)](@ref def_DensePolynomialZonotope)
* [Sparse polynomial zonotope (SparsePolynomialZonotope)](@ref def_SparsePolynomialZonotope)
* [Simplified sparse polynomial zonotope (SimpleSparsePolynomialZonotope)](@ref def_SimpleSparsePolynomialZonotope)
