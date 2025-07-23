```@contents
Pages = ["overapproximate_spz.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets.Approximations
```

# Overapproximation

```@docs
overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{SparsePolynomialZonotope}) where {N,S<:SparsePolynomialZonotope{N},
                                                                  NM,
                                                                  MAT<:MatrixZonotope{NM}}
overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{Zonotope}) where {N,S<:AbstractZonotope{N},NM,
                                                  MAT<:MatrixZonotope{NM}}
overapproximate(lm::LinearMap{N,S,NM,MAT},
                    ::Type{U}) where {N, S<:AbstractZonotope{N}, NM,
                                      MAT<:MatrixZonotopeProduct{NM},
                                      U<:Union{Zonotope, SparsePolynomialZonotope}}
_taylor_expmap
overapproximate(em::ExponentialMap{N,S,NM,MAT},
                                 ::Type{U},
                                 k::Int=2) where {N,
                                                  S<:Union{SparsePolynomialZonotope{N},
                                                           AbstractZonotope},
                                                  NM,
                                                  MAT<:AbstractMatrixZonotope{NM},
                                                  U<:Union{Zonotope,SparsePolynomialZonotope}}
```
