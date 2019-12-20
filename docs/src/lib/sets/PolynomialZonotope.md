```@meta
CurrentModule = LazySets
```

# [Polynomial zonotope (PolynomialZonotope)](@id def_PolynomialZonotope)

```@docs
PolynomialZonotope
dim(::PolynomialZonotope)
σ(::AbstractVector{N}, ::PolynomialZonotope{N}) where {N}
ρ(::AbstractVector{N}, ::PolynomialZonotope{N}) where {N}
polynomial_order(pz::PolynomialZonotope)
order(::PolynomialZonotope)
linear_map(::Matrix, ::PolynomialZonotope)
scale(::Number, ::PolynomialZonotope)
```
