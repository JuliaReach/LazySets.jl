```@meta
CurrentModule = LazySets
```

# [SparsePolynomialZonotope](@id def_SparsePolynomialZonotope)

```@docs
SparsePolynomialZonotope
rand(::Type{SparsePolynomialZonotope})
center(::SparsePolynomialZonotope)
genmat_dep(::SparsePolynomialZonotope)
genmat_indep(::SparsePolynomialZonotope)
expmat(::SparsePolynomialZonotope)
indexvector(P::SparsePolynomialZonotope)
uniqueID(::Int)
dim(::SparsePolynomialZonotope)
ngens_dep(::SparsePolynomialZonotope)
ngens_indep(::SparsePolynomialZonotope)
nparams(::SparsePolynomialZonotope)
order(::SparsePolynomialZonotope)
linear_map(::Union{Real, AbstractMatrix, LinearAlgebra.UniformScaling}, ::SparsePolynomialZonotope)
exact_sum(::SparsePolynomialZonotope, ::SparsePolynomialZonotope)
remove_redundant_generators(::SparsePolynomialZonotope)
```
