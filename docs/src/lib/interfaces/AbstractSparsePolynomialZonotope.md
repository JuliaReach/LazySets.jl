```@contents
Pages = ["AbstractSparsePolynomialZonotope.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Sparse polynomial zonotope sets (AbstractSparsePolynomialZonotope)](@id def_AbstractSparsePolynomialZonotope)

```@docs
AbstractSparsePolynomialZonotope
```

This interface requires to implement the following functions:

```@docs
expmat(::AbstractSparsePolynomialZonotope)
genmat_dep(::AbstractSparsePolynomialZonotope)
genmat_indep(::AbstractSparsePolynomialZonotope)
```

This interface defines the following functions:

```@docs
nparams(::AbstractSparsePolynomialZonotope)
```

## Implementations

* [Sparse polynomial zonotope (SparsePolynomialZonotope)](@ref def_SparsePolynomialZonotope)
* [Simplified sparse polynomial zonotope (SimpleSparsePolynomialZonotope)](@ref def_SimpleSparsePolynomialZonotope)
