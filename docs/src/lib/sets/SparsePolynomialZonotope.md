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
ngens_dep(::SparsePolynomialZonotope)
ngens_indep(::SparsePolynomialZonotope)
nparams(::SparsePolynomialZonotope)
order(::SparsePolynomialZonotope)
linear_map(::AbstractMatrix, ::SparsePolynomialZonotope)
translate(::SparsePolynomialZonotope, ::AbstractVector)
remove_redundant_generators(::SparsePolynomialZonotope)
reduce_order(::SparsePolynomialZonotope, ::Real, ::AbstractReductionMethod=GIR05())
ρ(::AbstractVector, ::SparsePolynomialZonotope)
```
Inherited from [`AbstractPolynomialZonotope`](@ref):
* [`dim`](@ref dim(::AbstractPolynomialZonotope))
