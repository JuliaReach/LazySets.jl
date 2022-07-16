```@meta
CurrentModule = LazySets
```

# [SimpleSparsePolynomialZonotope](@id def_SimpleSparsePolynomialZonotope)

```@docs
SimpleSparsePolynomialZonotope
PolynomialZonotope
center(::SimpleSparsePolynomialZonotope)
genmat(::SimpleSparsePolynomialZonotope)
expmat(::SimpleSparsePolynomialZonotope)
dim(::SimpleSparsePolynomialZonotope)
ngens(::SimpleSparsePolynomialZonotope)
nparams(::SimpleSparsePolynomialZonotope)
order(::SimpleSparsePolynomialZonotope)
linear_map(::AbstractMatrix, ::SimpleSparsePolynomialZonotope)
minkowski_sum(::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
cartesian_product(::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
linear_combination(::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
convex_hull(::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
convex_hull(::SimpleSparsePolynomialZonotope)
quadratic_map(::Vector{MT}, ::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
remove_redundant_generators(::SimpleSparsePolynomialZonotope)
```
