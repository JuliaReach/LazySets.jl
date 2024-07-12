```@meta
CurrentModule = LazySets
```

# [SparsePolynomialZonotope](@id def_SparsePolynomialZonotope)

```@docs
SparsePolynomialZonotope
center(::SparsePolynomialZonotope)
expmat(::SparsePolynomialZonotope)
genmat_dep(::SparsePolynomialZonotope)
genmat_indep(::SparsePolynomialZonotope)
indexvector(::SparsePolynomialZonotope)
ngens_dep(::SparsePolynomialZonotope)
ngens_indep(::SparsePolynomialZonotope)
polynomial_order(::SparsePolynomialZonotope)
rand(::Type{SparsePolynomialZonotope})
remove_redundant_generators(::SparsePolynomialZonotope)
uniqueID(::Int)
linear_map(::AbstractMatrix, ::SparsePolynomialZonotope)
reduce_order(::SparsePolynomialZonotope, ::Real, ::AbstractReductionMethod=GIR05())
ρ(::AbstractVector, ::SparsePolynomialZonotope)
translate(::SparsePolynomialZonotope, ::AbstractVector)
```

Inherited from [`AbstractPolynomialZonotope`](@ref):
* [`dim`](@ref dim(::AbstractPolynomialZonotope))
* [`nparams`](@ref dim(::AbstractPolynomialZonotope))
* [`order`](@ref dim(::AbstractPolynomialZonotope))
