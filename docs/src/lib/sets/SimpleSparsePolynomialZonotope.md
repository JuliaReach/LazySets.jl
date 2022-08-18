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
dim(::SimpleSparsePolynomialZonotope)
ngens(::SimpleSparsePolynomialZonotope)
nparams(::SimpleSparsePolynomialZonotope)
order(::SimpleSparsePolynomialZonotope)
linear_map(::Union{Real, AbstractMatrix, LinearAlgebra.UniformScaling}, ::SimpleSparsePolynomialZonotope)
minkowski_sum(::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
cartesian_product(::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
linear_combination(::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
convex_hull(::SimpleSparsePolynomialZonotope, ::SimpleSparsePolynomialZonotope)
convex_hull(::SimpleSparsePolynomialZonotope)
quadratic_map(::Vector{MT}, ::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
quadratic_map(Q::Vector{MT}, S1::SimpleSparsePolynomialZonotope, S2::SimpleSparsePolynomialZonotope) where {N, MT<:AbstractMatrix{N}}
remove_redundant_generators(::SimpleSparsePolynomialZonotope)
```
