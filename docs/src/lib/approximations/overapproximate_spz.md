```@contents
Pages = ["overapproximate_spz.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets.Approximations
```

# Overapproximation

```@docs
overapproximate(lm::LinearMap{N,S,NM,MAT}) where {N,S<:SparsePolynomialZonotope{N},NM,
                                                           MAT<:MatrixZonotope{NM}}
overapproximate(lm::LinearMap{N,S,NM,MAT}) where {N,S<:AbstractZonotope{N},NM,
                                                           MAT<:MatrixZonotope{NM}}
_taylor_expmap
overapproximate(::ExponentialMap, ::Int)
```
