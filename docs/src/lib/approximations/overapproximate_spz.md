```@contents
Pages = ["overapproximate_spz.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets.Approximations
```

# Overapproximation

```@docs
overapproximate(::LinearMap{N, SparsePolynomialZonotope{N}, NM, MAT}) where {N, NM,
	                 MAT <: MatrixZonotope{NM}}
_taylor_expmap
overapproximate(::ExponentialMap, ::Int)
```
