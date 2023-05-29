```@meta
CurrentModule = LazySets
```

# [SimpleSparsePolynomialZonotope](@id def_SimpleSparsePolynomialZonotope)

```@docs
SimpleSparsePolynomialZonotope
PolynomialZonotope
rand(::Type{SimpleSparsePolynomialZonotope})
center(::SimpleSparsePolynomialZonotope)
genmat(::SimpleSparsePolynomialZonotope)
expmat(::SimpleSparsePolynomialZonotope)
ngens(::SimpleSparsePolynomialZonotope)
nparams(::SimpleSparsePolynomialZonotope)
order(::SimpleSparsePolynomialZonotope)
linear_map(::AbstractMatrix, ::SimpleSparsePolynomialZonotope)
convex_hull(::SimpleSparsePolynomialZonotope)
quadratic_map(::Vector{MT}, ::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
quadratic_map(Q::Vector{MT}, S1::SimpleSparsePolynomialZonotope, S2::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
remove_redundant_generators(::SimpleSparsePolynomialZonotope)
```
Inherited from [`AbstractPolynomialZonotope`](@ref):
* [`dim`](@ref dim(::AbstractPolynomialZonotope))
