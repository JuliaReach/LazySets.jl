```@meta
CurrentModule = LazySets
```

# [SimpleSparsePolynomialZonotope](@id def_SimpleSparsePolynomialZonotope)

```@docs
SimpleSparsePolynomialZonotope
PolynomialZonotope
center(::SimpleSparsePolynomialZonotope)
convex_hull(::SimpleSparsePolynomialZonotope)
expmat(::SimpleSparsePolynomialZonotope)
genmat(::SimpleSparsePolynomialZonotope)
ngens(::SimpleSparsePolynomialZonotope)
nparams(::SimpleSparsePolynomialZonotope)
rand(::Type{SimpleSparsePolynomialZonotope})
remove_redundant_generators(::SimpleSparsePolynomialZonotope)
linear_map(::AbstractMatrix, ::SimpleSparsePolynomialZonotope)
quadratic_map(::Vector{MT}, ::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
quadratic_map(::Vector{MT}, ::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
```

Inherited from [`AbstractPolynomialZonotope`](@ref):
* [`dim`](@ref dim(::AbstractPolynomialZonotope))
* [`order`](@ref dim(::AbstractPolynomialZonotope))
